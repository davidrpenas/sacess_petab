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
namespace model_Isensee_JCB2018 {

void xdot_Isensee_JCB2018(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dACdt = -flux_v10_v_9 + flux_v12_v_11 - flux_v6_v_5 + flux_v7_v_6;  // xdot[0]
    dpACdt = flux_v10_v_9 - flux_v12_v_11 - flux_v8_v_7 + flux_v9_v_8;  // xdot[1]
    dAC_Fskdt = -flux_v11_v_10 + flux_v13_v_12 + flux_v6_v_5 - flux_v7_v_6;  // xdot[2]
    dpAC_Fskdt = flux_v11_v_10 - flux_v13_v_12 + flux_v8_v_7 - flux_v9_v_8;  // xdot[3]
    dcAMPdt = flux_v14_v_13 + flux_v15_v_14 - flux_v16_v_15;  // xdot[4]
    dCsubdt = flux_v17_v_16 - flux_v18_v_17 + flux_v19_v_18 - flux_v20_v_19 + flux_v27_v_26 + flux_v38_v_37 - flux_v40_v_39 + flux_v41_v_40;  // xdot[5]
    dCsub_H89dt = flux_v40_v_39 - flux_v41_v_40;  // xdot[6]
    dRII_C_2dt = -flux_v17_v_16 + flux_v18_v_17 - flux_v22_v_21 + flux_v23_v_22;  // xdot[7]
    dRIIp_C_2dt = -flux_v19_v_18 + flux_v20_v_19 + flux_v22_v_21 - flux_v23_v_22 - flux_v24_v_23 + flux_v25_v_24 + flux_v26_v_25 - flux_v30_v_29 + flux_v31_v_30 - flux_v32_v_31 + flux_v33_v_32 - flux_v34_v_33 + flux_v35_v_34 - flux_v36_v_35 + flux_v37_v_36;  // xdot[8]
    dRIIp_cAMP_C_2dt = flux_v24_v_23 - flux_v25_v_24 - flux_v26_v_25 - flux_v27_v_26;  // xdot[9]
    dRIIp_cAMP_2dt = flux_v27_v_26 - flux_v28_v_27 - flux_v29_v_28;  // xdot[10]
    dSp8_Br_cAMPSdt = flux_v1_v_0;  // xdot[11]
    dRp8_Br_cAMPSdt = flux_v2_v_1;  // xdot[12]
    dRp8_pCPT_cAMPSdt = flux_v3_v_2;  // xdot[13]
    dRp_cAMPSdt = flux_v4_v_3;  // xdot[14]
    dRIIp_Rp8_Br_cAMPS_C_2dt = flux_v30_v_29 - flux_v31_v_30;  // xdot[15]
    dRIIp_Rp8_pCPT_cAMPS_C_2dt = flux_v32_v_31 - flux_v33_v_32;  // xdot[16]
    dRIIp_Rp_cAMPS_C_2dt = flux_v34_v_33 - flux_v35_v_34;  // xdot[17]
    dRIIp_Sp8_Br_cAMPS_C_2dt = flux_v36_v_35 - flux_v37_v_36 - flux_v38_v_37;  // xdot[18]
    dRIIp_Sp8_Br_cAMPS_2dt = flux_v38_v_37 - flux_v39_v_38;  // xdot[19]
    dRIIp_2dt = flux_v19_v_18 - flux_v20_v_19 - flux_v21_v_20 + flux_v28_v_27 + flux_v29_v_28 + flux_v39_v_38;  // xdot[20]
    dRII_2dt = flux_v17_v_16 - flux_v18_v_17 + flux_v21_v_20;  // xdot[21]
    dIBMXdt = flux_v5_v_4;  // xdot[22]
    dPDEdt = -flux_v42_v_41 + flux_v43_v_42;  // xdot[23]
    dpPDEdt = flux_v42_v_41 - flux_v43_v_42;  // xdot[24]
}

} // namespace model_Isensee_JCB2018
} // namespace amici
