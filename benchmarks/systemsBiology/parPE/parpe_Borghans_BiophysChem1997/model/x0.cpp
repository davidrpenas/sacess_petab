#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void x0_Borghans_BiophysChem1997(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_Z_state;
    x0[1] = init_Y_state;
    x0[2] = init_A_state;
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
