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
namespace model_Bertozzi_PNAS2020 {

void sigmay_Bertozzi_PNAS2020(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_y_I_NY = sd_I_NY;  // sigmay[0]
    sigma_y_I_CA = sd_I_CA;  // sigmay[1]
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
