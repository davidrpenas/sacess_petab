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
namespace model_Brannmark_JBC2010 {

void xdot_Brannmark_JBC2010(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dIRdt = -flux_v1_v_0 + flux_v2_v_1 + flux_v6_v_5 + flux_v7_v_6;  // xdot[0]
    dIRinsdt = flux_v1_v_0 - flux_v2_v_1 - flux_v3_v_2;  // xdot[1]
    dIRpdt = flux_v3_v_2 - flux_v4_v_3 - flux_v6_v_5;  // xdot[2]
    dIRiPdt = flux_v4_v_3 - flux_v5_v_4;  // xdot[3]
    dIRidt = flux_v5_v_4 - flux_v7_v_6;  // xdot[4]
    dIRSdt = -flux_v8_v_7 + flux_v9_v_8;  // xdot[5]
    dIRSiPdt = flux_v8_v_7 - flux_v9_v_8;  // xdot[6]
    dXdt = -flux_v10_v_9 + flux_v11_v_10;  // xdot[7]
    dXpdt = flux_v10_v_9 - flux_v11_v_10;  // xdot[8]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
