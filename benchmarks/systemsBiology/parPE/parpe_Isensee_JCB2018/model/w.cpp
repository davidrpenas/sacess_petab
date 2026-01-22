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

namespace amici {
namespace model_Isensee_JCB2018 {

void w_Isensee_JCB2018(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    Fsk = Fsk_level*Heaviside_0;  // w[0]
    H89 = H89_level*Heaviside_2;  // w[1]
    IBMXex = Heaviside_4*IBMX_level;  // w[2]
    Rp8_Br_cAMPS_pAB = Heaviside_6*Rp8_Br_cAMPS_pAB_level;  // w[3]
    Rp8_pCPT_cAMPS_pAB = Heaviside_8*Rp8_pCPT_cAMPS_pAB_level;  // w[4]
    Rp_cAMPS_pAB = Heaviside_10*Rp_cAMPS_pAB_level;  // w[5]
    Sp8_Br_cAMPS_AM = Heaviside_12*Sp8_Br_cAMPS_AM_level;  // w[6]
    fourABnOH = Heaviside_14*fourABnOH_level;  // w[7]
    flux_v1_v_0 = -1.0*ki_Sp8_Br_cAMPS_AM*(Sp8_Br_cAMPS - Sp8_Br_cAMPS_AM*xi_i_Sp8_Br_cAMPS_AM);  // w[8]
    flux_v2_v_1 = -1.0*ki_Rp8_Br_cAMPS_pAB*(Rp8_Br_cAMPS - Rp8_Br_cAMPS_pAB*xi_i_Rp8_Br_cAMPS_pAB);  // w[9]
    flux_v3_v_2 = -1.0*ki_Rp8_pCPT_cAMPS_pAB*(Rp8_pCPT_cAMPS - Rp8_pCPT_cAMPS_pAB*xi_i_Rp8_pCPT_cAMPS_pAB);  // w[10]
    flux_v4_v_3 = -1.0*ki_Rp_cAMPS_pAB*(Rp_cAMPS - Rp_cAMPS_pAB*xi_i_Rp_cAMPS_pAB);  // w[11]
    flux_v5_v_4 = -1.0*ki_IBMX*(IBMX - IBMXex);  // w[12]
    flux_v6_v_5 = 1.0*AC*Fsk*kf_Fsk;  // w[13]
    flux_v7_v_6 = 1.0*AC_Fsk*KD_Fsk*kf_Fsk;  // w[14]
    flux_v8_v_7 = 1.0*Fsk*kf_Fsk*pAC;  // w[15]
    flux_v9_v_8 = 1.0*KD_Fsk*kf_Fsk*pAC_Fsk;  // w[16]
    flux_v10_v_9 = 1.0*AC*Csub*kp_AC;  // w[17]
    flux_v11_v_10 = 1.0*AC_Fsk*Csub*kp_AC;  // w[18]
    flux_v12_v_11 = 1.0*kdp_AC*pAC;  // w[19]
    flux_v13_v_12 = 1.0*kdp_AC*pAC_Fsk;  // w[20]
    flux_v14_v_13 = 1.0*ks_AC_cAMP*(AC + pAC*xi_pAC);  // w[21]
    flux_v15_v_14 = 1.0*ks_AC_cAMP*xi_AC_cAMP_Fsk*(AC_Fsk + pAC_Fsk*xi_pAC);  // w[22]
    flux_v16_v_15 = 1.0*KD_IBMX*cAMP*kdeg_cAMP_free*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // w[23]
    flux_v17_v_16 = 1.0*RII_C_2*kf_RII_C_2__RII_2;  // w[24]
    flux_v18_v_17 = 1.0*Csub*RII_2*kf_RII_2__RII_C_2;  // w[25]
    flux_v19_v_18 = 1.0*RIIp_C_2*kf_RII_C_2__RII_2*xi_kf_RII_C_2__RII_2;  // w[26]
    flux_v20_v_19 = 1.0*Csub*RIIp_2*kf_RII_2__RII_C_2*xi_kf_RII_2__RII_C_2;  // w[27]
    flux_v21_v_20 = 1.0*RIIp_2*kf_RIIp_2__RII_2;  // w[28]
    flux_v22_v_21 = 1.0*RII_C_2*kf_RII_C_2__RIIp_C_2;  // w[29]
    flux_v23_v_22 = 1.0*RIIp_C_2*kf_RIIp_C_2__RII_C_2;  // w[30]
    flux_v24_v_23 = 1.0*RIIp_C_2*cAMP*kf_cAMP;  // w[31]
    flux_v25_v_24 = 1.0*KD_cAMP*RIIp_cAMP_C_2*kf_cAMP;  // w[32]
    flux_v26_v_25 = 1.0*KD_IBMX*RIIp_cAMP_C_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // w[33]
    flux_v27_v_26 = 1.0*RIIp_cAMP_C_2*kf_RIIp_cAMP_C_2__RIIp_2;  // w[34]
    flux_v28_v_27 = 1.0*KD_cAMP*RIIp_cAMP_2*kf_cAMP;  // w[35]
    flux_v29_v_28 = 1.0*KD_IBMX*RIIp_cAMP_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // w[36]
    flux_v30_v_29 = 1.0*RIIp_C_2*Rp8_Br_cAMPS*kf_cAMP*xi_b_Rp8_Br_cAMPS;  // w[37]
    flux_v31_v_30 = 1.0*KD_cAMP*RIIp_Rp8_Br_cAMPS_C_2*kf_cAMP*xi_KD_Rp8_Br_cAMPS*xi_b_Rp8_Br_cAMPS;  // w[38]
    flux_v32_v_31 = 1.0*RIIp_C_2*Rp8_pCPT_cAMPS*kf_cAMP*xi_b_Rp8_pCPT_cAMPS;  // w[39]
    flux_v33_v_32 = 1.0*KD_cAMP*RIIp_Rp8_pCPT_cAMPS_C_2*kf_cAMP*xi_KD_Rp8_pCPT_cAMPS*xi_b_Rp8_pCPT_cAMPS;  // w[40]
    flux_v34_v_33 = 1.0*RIIp_C_2*Rp_cAMPS*kf_cAMP*xi_b_Rp_cAMPS;  // w[41]
    flux_v35_v_34 = 1.0*KD_cAMP*RIIp_Rp_cAMPS_C_2*kf_cAMP*xi_KD_Rp_cAMPS*xi_b_Rp_cAMPS;  // w[42]
    flux_v36_v_35 = 1.0*RIIp_C_2*Sp8_Br_cAMPS*kf_cAMP*xi_b_Sp8_Br_cAMPS;  // w[43]
    flux_v37_v_36 = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_C_2*kf_cAMP*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // w[44]
    flux_v38_v_37 = 1.0*RIIp_Sp8_Br_cAMPS_C_2*kf_RIIp_cAMP_C_2__RIIp_2;  // w[45]
    flux_v39_v_38 = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_2*kf_cAMP*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // w[46]
    flux_v40_v_39 = 1.0*Csub*H89*kf_H89;  // w[47]
    flux_v41_v_40 = 1.0*Csub_H89*KD_H89*kf_H89;  // w[48]
    flux_v42_v_41 = 1.0*Csub*PDE*kf_PDE_Csub;  // w[49]
    flux_v43_v_42 = 1.0*KD_PDE_Csub*kf_PDE_Csub*pPDE;  // w[50]
}

} // namespace model_Isensee_JCB2018
} // namespace amici
