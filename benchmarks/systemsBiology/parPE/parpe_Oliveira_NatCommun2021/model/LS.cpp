#include <parpeamici/optimizationApplication.h>
#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include <parpeoptimization/optimizationResultWriter.h>
#include <parpeoptimization/optimizer.h>
#include <parpecommon/logging.h>
#include <parpeoptimization/optimizationOptions.h>
#include <parpecommon/misc.h>
#include <chrono>
#include <amici/model.h>
#include <memory>
#include <string>
#include "cost_function.h"
//#include <mpi.h>
#include <cfloat>

using namespace std::chrono;
using namespace parpe;

namespace amici::generic_model {
    std::unique_ptr<amici::Model> getModel();
}


#ifdef __cplusplus
extern "C" {
#endif

class SilentLogger : public parpe::Logger {
public:
    SilentLogger() = default;

    void logmessage(parpe::loglevel lvl, std::string const& msg) const  {}
    void logmessage(parpe::loglevel lvl, const char *format, ...) const  {}
    void logmessage(parpe::loglevel lvl, const char *format, va_list argptr) const  {}
};

class CountingGradientFunction : public parpe::GradientFunction {
public:
    CountingGradientFunction(std::unique_ptr<parpe::GradientFunction> originalFunction, double init, double max)
        : originalFunction_(std::move(originalFunction)), numEvaluations_(0), init_time(init), max_time(max) {}

    FunctionEvaluationStatus evaluate(gsl::span<const double> parameters,
                                      double &fval,
                                      gsl::span<double> gradient,
                                      parpe::Logger *logger = nullptr,
                                      double *cpuTime = nullptr) const override {
        ++numEvaluations_; // Now atomic
/*        double current = MPI_Wtime();
        if ((current - init_time) >= max_time) {
            fval = std::numeric_limits<double>::max();
            std::fill(gradient.begin(), gradient.end(), std::numeric_limits<double>::max());
            return FunctionEvaluationStatus::functionEvaluationFailure;
        }
*/
        return originalFunction_->evaluate(parameters, fval, gradient, logger, cpuTime);
    }
	

    int getNumEvaluations() const { return numEvaluations_; }

    int numParameters() const override { return originalFunction_->numParameters(); }
    std::vector<std::string> getParameterIds() const override { return originalFunction_->getParameterIds(); }

private:
    std::unique_ptr<parpe::GradientFunction> originalFunction_;
    mutable std::atomic<int> numEvaluations_;
    double init_time; 
    double max_time; 
};

class MyOptimizationApplication : public parpe::OptimizationApplication {
public:
    using OptimizationApplication::OptimizationApplication;

    int id_local;
    int max_iters;
    std::vector<double> solution_vect;
    std::string outFileArgument_save;
    std::string outpath;
    int numEvaluations;
    double init_time;
    double max_time;
	
    OptimizationProblem* getProblem() {
        return problem.get();
    }

    virtual void finalizeTiming(double wallTimeSeconds, double cpuTimeSeconds) override {}

	virtual void initProblem(const std::string& inFileArgument, const std::string& outFileArgument) override {
		outFileArgument_save = outpath + "/temp_file_parPE_" + std::to_string(id_local) + ".h5";

		auto h5Outfile = parpe::hdf5CreateFile(outFileArgument_save, true);
		logParPEVersion(h5Outfile);

		dataProvider = std::make_unique<parpe::MultiConditionDataProviderHDF5>(
			amici::generic_model::getModel(), inFileArgument);

		auto h5Infile = dataProvider->getHdf5File();
		auto optimizationOptions = parpe::OptimizationOptions::fromHDF5(h5Infile);



		auto multiCondProb = new parpe::MultiConditionProblem(
		     dataProvider.get(), 
	  	      nullptr, 
		     std::make_unique<SilentLogger>(), 
	             std::make_unique<parpe::OptimizationResultWriter>(h5Outfile, "/multistarts/")
		);

		auto originalFunction = multiCondProb->cost_fun_.release();
		auto countingFunction = std::make_unique<CountingGradientFunction>(
			std::unique_ptr<parpe::GradientFunction>(originalFunction), init_time, max_time
		);
		multiCondProb->cost_fun_ = std::move(countingFunction);

		optimizationOptions->optimizer = parpe::optimizerName::OPTIMIZER_IPOPT;
		optimizationOptions->maxOptimizerIterations = max_iters;
		optimizationOptions->printToStdout = false;
		optimizationOptions->multistartsInParallel = false;
		optimizationOptions->numStarts = 0;

		multiCondProb->setInitialParameters(solution_vect);
		multiCondProb->setOptimizationOptions(*optimizationOptions);

		problem.reset(multiCondProb);
	}


	int getNumEvaluations() const {
		auto *countingFunction = dynamic_cast<CountingGradientFunction*>(problem->cost_fun_.get());
		if (countingFunction) {
			return countingFunction->getNumEvaluations();
		}
		return -1; 
	}
	

    virtual ~MyOptimizationApplication() override {
        remove( outFileArgument_save.c_str() );
		withMPI = false;
    }

    virtual void printUsage(char* const argZero) override {}

    virtual void logParPEVersion(const H5::H5File &file) const override {}
	
private:
    /** DataProvider as interface to HDF5 data */
    std::unique_ptr<parpe::MultiConditionDataProviderHDF5> dataProvider;
};

int call_parPE(int id, double *fx, int N, int *nevals, int iters, double init_time, double max_time, double *sol_vect_ptr, const char *name_input, const char *name_output) {

#ifndef NDEBUG
    // Set stdout to unbuffered when debugging
    setbuf(stdout, NULL);
#endif
    setenv("PARPE_LOG_SIMULATIONS","0",1);
    setenv("PARPE_NO_DEBUG","1",1);
    parpe::minimumLogLevel = loglevel::critical;

    std::vector<double> sol_vect(N);
    for (int i=0;i<N;++i) sol_vect[i] = sol_vect_ptr[i];

    MyOptimizationApplication app;
	
//    const char *files[3];
//    files[0] = const_cast<char*>("program");
//    files[1] = name_input;
//    files[2] = name_output;

    app.id_local = id;
    app.solution_vect = sol_vect;
    app.outpath = name_output;
    app.max_iters = iters;
    app.init_time = init_time;
    app.max_time = max_time;


    char ** files = (char **) malloc(3*sizeof(char *));
    files[0] = "program";
    files[1] = (char *) name_input;
    files[2] = "output";

    app.id_local = id;
    app.solution_vect = sol_vect;
    app.outpath = name_output;

    std::tuple<int, double, std::vector<double, std::allocator<double>>> result = app.run(2, files);
    *nevals = *nevals + app.getNumEvaluations();

    int  status = std::get<0>(result);
    *fx = std::get<1>(result);
    std::vector<double> end_solution = std::get<2>(result);

    
    for (int i=0;i<N;++i) sol_vect_ptr[i] = end_solution[i];

   
    free(files);
    end_solution.clear();
    sol_vect.clear();

    return 1;
}

#ifdef __cplusplus
}
#endif
