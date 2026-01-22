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
namespace model_Giordano_Nature2020 {

void xdot_Giordano_Nature2020(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dSusceptibledt = -flux_Susceptible_to_Infected;  // xdot[0]
    dInfecteddt = -flux_Infected_to_Ailing - flux_Infected_to_Diagnosed - flux_Infected_to_Healed + flux_Susceptible_to_Infected;  // xdot[1]
    dDiagnoseddt = -flux_Diagnosed_to_Healed - flux_Diagnosed_to_Recognized + flux_Infected_to_Diagnosed;  // xdot[2]
    dAilingdt = -flux_Ailing_to_Healed - flux_Ailing_to_Recognised - flux_Ailing_to_Threatened + flux_Infected_to_Ailing;  // xdot[3]
    dRecognizeddt = flux_Ailing_to_Recognised + flux_Diagnosed_to_Recognized - flux_Recognised_to_Healed - flux_Recognised_to_Threatened;  // xdot[4]
    dThreateneddt = flux_Ailing_to_Threatened + flux_Recognised_to_Threatened - flux_Threatened_to_Extinct - flux_Threatened_to_Healed;  // xdot[5]
    dHealeddt = flux_Ailing_to_Healed + flux_Diagnosed_to_Healed + flux_Infected_to_Healed + flux_Recognised_to_Healed + flux_Threatened_to_Healed;  // xdot[6]
    dExtinctdt = flux_Threatened_to_Extinct;  // xdot[7]
    dCumulativeDiagnoseddt = flux_Ailing_to_Recognised + flux_Ailing_to_Threatened + flux_Infected_to_Diagnosed;  // xdot[8]
    dDiagnosedHealeddt = flux_Diagnosed_to_Healed + flux_Recognised_to_Healed + flux_Threatened_to_Healed;  // xdot[9]
}

} // namespace model_Giordano_Nature2020
} // namespace amici
