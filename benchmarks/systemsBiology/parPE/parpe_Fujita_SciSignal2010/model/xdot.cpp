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
namespace model_Fujita_SciSignal2010 {

void xdot_Fujita_SciSignal2010(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dEGFRdt = -flux_v10_reaction_10 + flux_v11_reaction_11 - flux_v1_reaction_1;  // xdot[0]
    dpEGFRdt = -flux_v2_reaction_2 + flux_v3_reaction_3 - flux_v4_reaction_4 + flux_v9_reaction_9;  // xdot[1]
    dpEGFR_Aktdt = flux_v2_reaction_2 - flux_v3_reaction_3;  // xdot[2]
    dAktdt = -flux_v2_reaction_2 + flux_v7_reaction_7;  // xdot[3]
    dpAktdt = flux_v3_reaction_3 - flux_v5_reaction_5 + flux_v6_reaction_6 - flux_v7_reaction_7;  // xdot[4]
    dS6dt = -flux_v5_reaction_5 + flux_v8_reaction_8;  // xdot[5]
    dpAkt_S6dt = flux_v5_reaction_5 - flux_v6_reaction_6;  // xdot[6]
    dpS6dt = flux_v6_reaction_6 - flux_v8_reaction_8;  // xdot[7]
    dEGF_EGFRdt = flux_v1_reaction_1 - flux_v9_reaction_9;  // xdot[8]
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
