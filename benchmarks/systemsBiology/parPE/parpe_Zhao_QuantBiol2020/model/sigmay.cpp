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
namespace model_Zhao_QuantBiol2020 {

void sigmay_Zhao_QuantBiol2020(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_confirmed_infected = noiseParameter1_observable_confirmed_infected;  // sigmay[0]
    sigma_observable_susceptible = noiseParameter1_observable_susceptible;  // sigmay[1]
    sigma_observable_unquarantined_infected = noiseParameter1_observable_confirmed_infected;  // sigmay[2]
    sigma_observable_quarantined_infected = noiseParameter1_observable_quarantined_infected;  // sigmay[3]
    sigma_observable_cumulative_infected = noiseParameter1_observable_cumulative_infected;  // sigmay[4]
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
