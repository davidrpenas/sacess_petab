#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Rahman_MBS2016 {

void x_solver_Rahman_MBS2016(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = susceptible;
    x_solver[1] = infected_normal;
    x_solver[2] = infected_moderate;
    x_solver[3] = infected_weak;
    x_solver[4] = treated_normal;
    x_solver[5] = treated_moderate;
    x_solver[6] = treated_weak;
}

} // namespace model_Rahman_MBS2016
} // namespace amici
