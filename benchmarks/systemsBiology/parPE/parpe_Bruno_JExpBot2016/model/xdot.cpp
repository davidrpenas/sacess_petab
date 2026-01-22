#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void xdot_Bruno_JExpBot2016(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dbcardt = -flux_v1_ReactionName;  // xdot[0]
    dbcrydt = -flux_v3_ReactionName - flux_v4_ReactionName;  // xdot[1]
    db10dt = flux_v1_ReactionName - flux_v2_ReactionName + flux_v3_ReactionName;  // xdot[2]
    dbiodt = flux_v1_ReactionName + flux_v2_ReactionName + flux_v4_ReactionName;  // xdot[3]
    dohb10dt = flux_v4_ReactionName - flux_v5_ReactionName + flux_v6_ReactionName;  // xdot[4]
    dohbiodt = flux_v3_ReactionName + flux_v5_ReactionName + flux_v6_ReactionName;  // xdot[5]
    dzeadt = -flux_v6_ReactionName;  // xdot[6]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
