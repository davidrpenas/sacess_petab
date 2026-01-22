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
namespace model_Weber_BMC2015 {

void xdot_Weber_BMC2015(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dPKDdt = -flux_v1_v_0 - flux_v2_v_1 + flux_v3_v_2 + flux_v4_v_3 + flux_v5_v_4 - flux_v6_v_5;  // xdot[0]
    dPKDDAGadt = flux_v1_v_0 + flux_v2_v_1 - flux_v3_v_2 - flux_v7_v_6;  // xdot[1]
    dPI4K3Bdt = flux_v10_v_9 + flux_v11_v_10 - flux_v12_v_11 + flux_v8_v_7 - flux_v9_v_8;  // xdot[2]
    dPI4K3Badt = -flux_v13_v_12 - flux_v8_v_7 + flux_v9_v_8;  // xdot[3]
    dCERTERadt = -flux_v14_v_13 + flux_v15_v_14 + flux_v17_v_16 + flux_v18_v_17 - flux_v19_v_18;  // xdot[4]
    dCERTdt = -flux_v15_v_14 + flux_v16_v_15 - flux_v20_v_19;  // xdot[5]
    dCERTTGNadt = flux_v14_v_13 - flux_v16_v_15 - flux_v21_v_20;  // xdot[6]
}

} // namespace model_Weber_BMC2015
} // namespace amici
