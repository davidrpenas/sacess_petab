#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Sneyd_PNAS2002 {

void x_solver_Sneyd_PNAS2002(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = IPR_O;
    x_solver[1] = IPR_R;
    x_solver[2] = IPR_I1;
    x_solver[3] = IPR_S;
    x_solver[4] = IPR_A;
    x_solver[5] = IPR_I2;
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
