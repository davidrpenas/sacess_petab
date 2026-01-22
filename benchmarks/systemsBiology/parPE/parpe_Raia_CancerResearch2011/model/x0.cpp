#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void x0_Raia_CancerResearch2011(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1.3;
    x0[1] = init_Rec_i;
    x0[5] = 2.7999999999999998;
    x0[7] = 165.0;
    x0[10] = 0.34000000000000002;
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
