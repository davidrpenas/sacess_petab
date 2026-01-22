#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void x0_Isensee_JCB2018(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1.0;
    x0[4] = 0.057521875897794901;
    x0[5] = 0.19205140543353799;
    x0[7] = 0.51456221322303897;
    x0[8] = 0.28797435264320298;
    x0[9] = 0.0054120287002202899;
    x0[10] = 0.00040041663400605898;
    x0[20] = 0.13397950594491601;
    x0[21] = 0.057671482854615999;
    x0[23] = 1.0;
}

} // namespace model_Isensee_JCB2018
} // namespace amici
