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
namespace model_Elowitz_Nature2000 {

void xdot_Elowitz_Nature2000(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dX_proteindt = flux_v5_v_4 - flux_v9_v_8;  // xdot[0]
    dY_proteindt = -flux_v10_v_9 + flux_v6_v_5;  // xdot[1]
    dZ_proteindt = -flux_v11_v_10 + flux_v7_v_6;  // xdot[2]
    dGFPdt = -flux_v12_v_11 + flux_v16_v_15 + flux_v8_v_7;  // xdot[3]
    dX_mRNAdt = flux_v13_v_12 - flux_v1_v_0;  // xdot[4]
    dY_mRNAdt = flux_v14_v_13 - flux_v2_v_1;  // xdot[5]
    dZ_mRNAdt = flux_v15_v_14 - flux_v3_v_2;  // xdot[6]
    dGFP_mRNAdt = -flux_v4_v_3;  // xdot[7]
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
