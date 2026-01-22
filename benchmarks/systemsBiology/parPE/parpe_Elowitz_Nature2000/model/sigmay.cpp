#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void sigmay_Elowitz_Nature2000(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_fluorescence = noiseParameter1_fluorescence;  // sigmay[0]
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
