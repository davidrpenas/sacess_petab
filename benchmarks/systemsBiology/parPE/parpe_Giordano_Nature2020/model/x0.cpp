#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void x0_Giordano_Nature2020(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.99999629999999995;
    x0[1] = 3.3333333299999999e-6;
    x0[2] = 3.33333333e-7;
    x0[3] = 1.6666666599999999e-8;
    x0[4] = 3.33333333e-8;
    x0[8] = 3.6666666630000003e-7;
}

} // namespace model_Giordano_Nature2020
} // namespace amici
