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
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void sigmay_Okuonghae_ChaosSolitonsFractals2020(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_Cumulative = noiseParameter1_observable_Cumulative;  // sigmay[0]
    sigma_observable_ActiveCases = noiseParameter1_observable_ActiveCases;  // sigmay[1]
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
