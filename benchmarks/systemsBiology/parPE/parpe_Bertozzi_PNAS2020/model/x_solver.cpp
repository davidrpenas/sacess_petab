#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void x_solver_Bertozzi_PNAS2020(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = I_;
    x_solver[1] = R_;
    x_solver[2] = S_;
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
