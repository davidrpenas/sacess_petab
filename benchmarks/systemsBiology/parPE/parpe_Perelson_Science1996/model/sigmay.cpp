#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Perelson_Science1996 {

void sigmay_Perelson_Science1996(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_task0_model0_perelson1_V = noiseParameter1_task0_model0_perelson1_V;  // sigmay[0]
}

} // namespace model_Perelson_Science1996
} // namespace amici
