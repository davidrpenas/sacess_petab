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
namespace model_Giordano_Nature2020 {

void sigmay_Giordano_Nature2020(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_TotalCases = noiseParameter1_observable_TotalCases;  // sigmay[0]
    sigma_observable_Deaths = noiseParameter1_observable_Deaths;  // sigmay[1]
    sigma_observable_Recovered = noiseParameter1_observable_Recovered;  // sigmay[2]
    sigma_observable_CurrentCases = noiseParameter1_observable_CurrentCases;  // sigmay[3]
    sigma_observable_DiagHome = noiseParameter1_observable_DiagHome;  // sigmay[4]
    sigma_observable_Hospit = noiseParameter1_observable_Hospit;  // sigmay[5]
    sigma_observable_ICU = noiseParameter1_observable_ICU;  // sigmay[6]
}

} // namespace model_Giordano_Nature2020
} // namespace amici
