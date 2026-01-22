#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void xdot_Oliveira_NatCommun2021(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dExposeddt = -flux_Exposed_to_Asymptomatic - flux_Exposed_to_Symptomatic + flux_Susceptible_to_Exposed;  // xdot[0]
    dSusceptibledt = -flux_Susceptible_to_Exposed;  // xdot[1]
    dAsymptomaticdt = -flux_Asymptomatic_to_Recovered + flux_Exposed_to_Asymptomatic;  // xdot[2]
    dSymptomaticdt = flux_Exposed_to_Symptomatic - flux_Symptomatic_to_Hospitalized - flux_Symptomatic_to_ICU - flux_Symptomatic_to_Recovered;  // xdot[3]
    dRecovereddt = flux_Asymptomatic_to_Recovered + flux_Hospitalized_to_Recovered + flux_Symptomatic_to_Recovered;  // xdot[4]
    dHospitaldt = -flux_Hospitalized_to_Deaths - flux_Hospitalized_to_ICU - flux_Hospitalized_to_Recovered + flux_ICU_to_Hospitalized + flux_Symptomatic_to_Hospitalized;  // xdot[5]
    dICUdt = flux_Hospitalized_to_ICU - flux_ICU_to_Deaths - flux_ICU_to_Hospitalized + flux_Symptomatic_to_ICU;  // xdot[6]
    dDeathsdt = flux_Hospitalized_to_Deaths + flux_ICU_to_Deaths;  // xdot[7]
    dCumulative_casesdt = 1.0*Exposed*exposed_to_symptomatic;  // xdot[8]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
