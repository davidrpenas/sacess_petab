#include <parpeamici/multiConditionDataProvider.h>
#include <parpeamici/multiConditionProblem.h>
#include <amici/model.h>
#include <amici/hdf5.h>
#include <iostream>
#include <memory>
#include <parpecommon/logging.h>

namespace amici::generic_model {
        extern std::unique_ptr<amici::Model> getModel();
}

using namespace parpe;

#ifdef __cplusplus
extern "C" {
#endif

// GLOBAL VARIABLES
std::vector<double> lb_Problem;
std::vector<double> ub_Problem;
int size_problem_Problem;
std::vector<double> optimizationParams_Problem;

MultiConditionProblem  *problem_g_Problem;
MultiConditionDataProviderHDF5 *dataProvider_g_Problem;

// IMPL...
void load_proplem(const char* input_string) {

    std::string inFileArgument = input_string;
    //parpe::logmessage(parpe::LOGLVL_INFO, "Reading options and data from '%s'.",  inFileArgument.c_str());
    std::string outFileArgument = "";

    // setup data and problem
    static MultiConditionDataProviderHDF5 dataProvider_Problem(amici::generic_model::getModel(), inFileArgument);

    dataProvider_g_Problem = &dataProvider_Problem;


    auto options = OptimizationOptions::fromHDF5(dataProvider_Problem.getHdf5File());


    std::unique_ptr<OptimizationResultWriter> rw;
    if (!outFileArgument.empty()) {
        rw = std::make_unique<OptimizationResultWriter>(outFileArgument, true, "/");
    }

    static MultiConditionProblem  problem( &dataProvider_Problem, nullptr, nullptr, std::move(rw));
    problem_g_Problem = &problem;

    // problem = problem_aux;

    // Read nominal parameters
    optimizationParams_Problem = amici::hdf5::getDoubleDataset1D(dataProvider_Problem.getHdf5File(), "/parameters/nominalValues");


    std::vector<double> lb(dataProvider_Problem.getNumOptimizationParameters());
    std::vector<double> ub(dataProvider_Problem.getNumOptimizationParameters());

    dataProvider_Problem.getOptimizationParametersUpperBounds(ub);
    dataProvider_Problem.getOptimizationParametersLowerBounds(lb);


    // fill global variables
    lb_Problem = lb;
    ub_Problem = ub;
    size_problem_Problem = dataProvider_Problem.getNumOptimizationParameters();

    double fval = NAN;
    for (int i=0; i < size_problem_Problem; ++i) optimizationParams_Problem[i] = 0;
    std::vector<double> gradient(optimizationParams_Problem.size(), NAN);
    std::vector<double> res;
    problem.cost_fun_->evaluate(optimizationParams_Problem, fval, gradient);

}


void return_bounds(double* ub, double* lb) {
    for (int i = 0; i < size_problem_Problem; ++i) {
        ub[i] = ub_Problem[i];
        lb[i] = lb_Problem[i];
    }
}

int return_size_problem() {
    return size_problem_Problem;
}



void eval_cost(double* solution, double* fx, double **res, int *size_res, double* grad) {

//    std::vector<double> optimizationParams_Problem(size_problem_Problem);
//    for (int i = 0; i < size_problem_Problem; ++i) {
//        optimizationParams_Problem[i] = solution[i];
//    }

     setenv("PARPE_LOG_SIMULATIONS","0",1);
     setenv("PARPE_NO_DEBUG","1",1);
     parpe::minimumLogLevel = parpe::loglevel::critical;

    double fval = NAN;
    std::vector<double> gradient(size_problem_Problem, NAN);
    std::vector<double> resi;

    
    for (int i=0; i < size_problem_Problem; ++i) optimizationParams_Problem[i] = solution[i]; 

//    problem_g_Problem->cost_fun_->evaluate(optimizationParams_Problem, fval, gradient);
    problem_g_Problem->cost_fun_->evaluate(optimizationParams_Problem, fval, gsl::span<double>());

    //results
    *fx = fval;

    for (int i = 0; i < size_problem_Problem; ++i) {
        grad[i] = gradient[i];
    }

	
    if (!resi.empty()) {
        *size_res = static_cast<int>(resi.size());
        *res = (double*)malloc(*size_res * sizeof(double));
        if (*res == nullptr) {
            perror("Memory allocation failed for residuals");
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < *size_res; ++i) {
            (*res)[i] = resi[i];
        }
    } else {
        *size_res = 0;
        *res = nullptr;
    }

}

#ifdef __cplusplus
}
#endif

//}
//}

