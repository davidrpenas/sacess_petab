#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void xdot_Zhao_QuantBiol2020(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dSusceptibledt = -flux_Susceptible_to_Unquarantined;  // xdot[0]
    dUnquarantined_Infecteddt = flux_Susceptible_to_Unquarantined - flux_Unquarantined_to_Quarantined;  // xdot[1]
    dQuarantined_Infecteddt = -flux_Quarantined_to_Confirmed + flux_Unquarantined_to_Quarantined;  // xdot[2]
    dConfirmed_Infecteddt = flux_Quarantined_to_Confirmed;  // xdot[3]
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
