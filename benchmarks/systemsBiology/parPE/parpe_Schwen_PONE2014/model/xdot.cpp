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
namespace model_Schwen_PONE2014 {

void xdot_Schwen_PONE2014(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dInsdt = -flux_v1_ka1 - flux_v2_ka2fold - flux_v3_v_2 + flux_v4_v_3 + flux_v5_kd1 + flux_v6_kd2fold;  // xdot[0]
    dRec1dt = flux_v11_v_10 - flux_v1_ka1 + flux_v5_kd1;  // xdot[1]
    dRec2dt = flux_v12_v_11 - flux_v2_ka2fold + flux_v6_kd2fold;  // xdot[2]
    dIR1dt = flux_v1_ka1 - flux_v5_kd1 - flux_v7_v_6 + flux_v9_v_8;  // xdot[3]
    dIR2dt = flux_v10_v_9 + flux_v2_ka2fold - flux_v6_kd2fold - flux_v8_v_7;  // xdot[4]
    dIR1indt = -flux_v11_v_10 + flux_v7_v_6 - flux_v9_v_8;  // xdot[5]
    dIR2indt = -flux_v10_v_9 - flux_v12_v_11 + flux_v8_v_7;  // xdot[6]
    dUptake1dt = flux_v13_v_12;  // xdot[7]
    dUptake2dt = flux_v14_v_13;  // xdot[8]
    dInsulinFragmentsdt = flux_v11_v_10 + flux_v12_v_11;  // xdot[9]
    dBoundUnspecdt = flux_v3_v_2 - flux_v4_v_3;  // xdot[10]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
