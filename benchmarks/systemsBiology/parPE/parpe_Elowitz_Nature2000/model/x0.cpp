#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void x0_Elowitz_Nature2000(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_X_protein;
    x0[1] = init_Y_protein;
    x0[2] = init_Z_protein;
    x0[3] = init_GFP;
    x0[4] = init_X_mRNA;
    x0[5] = init_Y_mRNA;
    x0[6] = init_Z_mRNA;
    x0[7] = init_GFP_mRNA;
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
