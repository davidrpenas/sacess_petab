#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void x_solver_Borghans_BiophysChem1997(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Z_state;
    x_solver[1] = Y_state;
    x_solver[2] = A_state;
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
