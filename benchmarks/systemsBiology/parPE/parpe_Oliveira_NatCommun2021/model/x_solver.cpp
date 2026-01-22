#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void x_solver_Oliveira_NatCommun2021(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Exposed;
    x_solver[1] = Susceptible;
    x_solver[2] = Asymptomatic;
    x_solver[3] = Symptomatic;
    x_solver[4] = Recovered;
    x_solver[5] = Hospital;
    x_solver[6] = ICU;
    x_solver[7] = Deaths;
    x_solver[8] = Cumulative_cases;
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
