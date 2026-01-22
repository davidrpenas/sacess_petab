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
namespace model_Oliveira_NatCommun2021 {

void sigmay_Oliveira_NatCommun2021(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_cumulative_deaths = noiseParameter1_cumulative_deaths;  // sigmay[0]
    sigma_cumulative_cases = noiseParameter1_cumulative_cases;  // sigmay[1]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
