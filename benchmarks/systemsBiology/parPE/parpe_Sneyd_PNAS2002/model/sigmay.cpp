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
namespace model_Sneyd_PNAS2002 {

void sigmay_Sneyd_PNAS2002(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_open_probability = noiseParameter1_open_probability;  // sigmay[0]
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
