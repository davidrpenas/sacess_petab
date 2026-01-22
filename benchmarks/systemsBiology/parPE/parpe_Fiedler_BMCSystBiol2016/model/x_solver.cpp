#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void x_solver_Fiedler_BMCSystBiol2016(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = RAF;
    x_solver[1] = pRAF;
    x_solver[2] = MEK;
    x_solver[3] = pMEK;
    x_solver[4] = ERK;
    x_solver[5] = pERK;
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
