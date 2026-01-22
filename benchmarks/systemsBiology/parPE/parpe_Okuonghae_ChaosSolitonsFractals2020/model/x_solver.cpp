#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void x_solver_Okuonghae_ChaosSolitonsFractals2020(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = susceptible;
    x_solver[1] = exposed;
    x_solver[2] = asymptomatic;
    x_solver[3] = symptomatic;
    x_solver[4] = detected;
    x_solver[5] = recovered;
    x_solver[6] = deceased;
    x_solver[7] = detected_cumulative;
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
