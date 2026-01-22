#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 26> dwdx_colptrs_Isensee_JCB2018_ = {
    0, 3, 6, 9, 12, 14, 20, 21, 23, 30, 33, 35, 37, 39, 41, 43, 44, 45, 46, 48, 49, 51, 52, 56, 60, 64
};

void dwdx_colptrs_Isensee_JCB2018(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Isensee_JCB2018_));
}
} // namespace model_Isensee_JCB2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 64> dwdx_rowvals_Isensee_JCB2018_ = {
    13, 17, 21, 15, 19, 21, 14, 18, 22, 16, 20, 22, 23, 31, 17, 18, 25, 27, 47, 49, 48, 24, 29, 26, 30, 31, 37, 39, 41, 43, 32, 33, 34, 35, 36, 8, 43, 9, 37, 10, 39, 11, 41, 38, 40, 42, 44, 45, 46, 27, 28, 25, 12, 23, 33, 36, 23, 33, 36, 49, 23, 33, 36, 50
};

void dwdx_rowvals_Isensee_JCB2018(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Isensee_JCB2018_));
}
} // namespace model_Isensee_JCB2018
} // namespace amici




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
#include "dwdx.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void dwdx_Isensee_JCB2018(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v14_v_13_dAC = 1.0*ks_AC_cAMP;  // dwdx[2]
        dflux_v12_v_11_dpAC = 1.0*kdp_AC;  // dwdx[4]
        dflux_v14_v_13_dpAC = 1.0*ks_AC_cAMP*xi_pAC;  // dwdx[5]
        dflux_v7_v_6_dAC_Fsk = 1.0*KD_Fsk*kf_Fsk;  // dwdx[6]
        dflux_v15_v_14_dAC_Fsk = 1.0*ks_AC_cAMP*xi_AC_cAMP_Fsk;  // dwdx[8]
        dflux_v9_v_8_dpAC_Fsk = 1.0*KD_Fsk*kf_Fsk;  // dwdx[9]
        dflux_v13_v_12_dpAC_Fsk = 1.0*kdp_AC;  // dwdx[10]
        dflux_v15_v_14_dpAC_Fsk = 1.0*ks_AC_cAMP*xi_AC_cAMP_Fsk*xi_pAC;  // dwdx[11]
        dflux_v41_v_40_dCsub_H89 = 1.0*KD_H89*kf_H89;  // dwdx[20]
        dflux_v17_v_16_dRII_C_2 = 1.0*kf_RII_C_2__RII_2;  // dwdx[21]
        dflux_v22_v_21_dRII_C_2 = 1.0*kf_RII_C_2__RIIp_C_2;  // dwdx[22]
        dflux_v19_v_18_dRIIp_C_2 = 1.0*kf_RII_C_2__RII_2*xi_kf_RII_C_2__RII_2;  // dwdx[23]
        dflux_v23_v_22_dRIIp_C_2 = 1.0*kf_RIIp_C_2__RII_C_2;  // dwdx[24]
        dflux_v25_v_24_dRIIp_cAMP_C_2 = 1.0*KD_cAMP*kf_cAMP;  // dwdx[30]
        dflux_v27_v_26_dRIIp_cAMP_C_2 = 1.0*kf_RIIp_cAMP_C_2__RIIp_2;  // dwdx[32]
        dflux_v28_v_27_dRIIp_cAMP_2 = 1.0*KD_cAMP*kf_cAMP;  // dwdx[33]
        dflux_v1_v_0_dSp8_Br_cAMPS = -1.0*ki_Sp8_Br_cAMPS_AM;  // dwdx[35]
        dflux_v2_v_1_dRp8_Br_cAMPS = -1.0*ki_Rp8_Br_cAMPS_pAB;  // dwdx[37]
        dflux_v3_v_2_dRp8_pCPT_cAMPS = -1.0*ki_Rp8_pCPT_cAMPS_pAB;  // dwdx[39]
        dflux_v4_v_3_dRp_cAMPS = -1.0*ki_Rp_cAMPS_pAB;  // dwdx[41]
        dflux_v31_v_30_dRIIp_Rp8_Br_cAMPS_C_2 = 1.0*KD_cAMP*kf_cAMP*xi_KD_Rp8_Br_cAMPS*xi_b_Rp8_Br_cAMPS;  // dwdx[43]
        dflux_v33_v_32_dRIIp_Rp8_pCPT_cAMPS_C_2 = 1.0*KD_cAMP*kf_cAMP*xi_KD_Rp8_pCPT_cAMPS*xi_b_Rp8_pCPT_cAMPS;  // dwdx[44]
        dflux_v35_v_34_dRIIp_Rp_cAMPS_C_2 = 1.0*KD_cAMP*kf_cAMP*xi_KD_Rp_cAMPS*xi_b_Rp_cAMPS;  // dwdx[45]
        dflux_v37_v_36_dRIIp_Sp8_Br_cAMPS_C_2 = 1.0*KD_cAMP*kf_cAMP*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // dwdx[46]
        dflux_v38_v_37_dRIIp_Sp8_Br_cAMPS_C_2 = 1.0*kf_RIIp_cAMP_C_2__RIIp_2;  // dwdx[47]
        dflux_v39_v_38_dRIIp_Sp8_Br_cAMPS_2 = 1.0*KD_cAMP*kf_cAMP*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // dwdx[48]
        dflux_v21_v_20_dRIIp_2 = 1.0*kf_RIIp_2__RII_2;  // dwdx[50]
        dflux_v5_v_4_dIBMX = -1.0*ki_IBMX;  // dwdx[52]
        dflux_v43_v_42_dpPDE = 1.0*KD_PDE_Csub*kf_PDE_Csub;  // dwdx[63]
    }

    // dynamic expressions
    dflux_v6_v_5_dAC = 1.0*Fsk*kf_Fsk;  // dwdx[0]
    dflux_v10_v_9_dAC = 1.0*Csub*kp_AC;  // dwdx[1]
    dflux_v8_v_7_dpAC = 1.0*Fsk*kf_Fsk;  // dwdx[3]
    dflux_v11_v_10_dAC_Fsk = 1.0*Csub*kp_AC;  // dwdx[7]
    dflux_v16_v_15_dcAMP = 1.0*KD_IBMX*kdeg_cAMP_free*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdx[12]
    dflux_v24_v_23_dcAMP = 1.0*RIIp_C_2*kf_cAMP;  // dwdx[13]
    dflux_v10_v_9_dCsub = 1.0*AC*kp_AC;  // dwdx[14]
    dflux_v11_v_10_dCsub = 1.0*AC_Fsk*kp_AC;  // dwdx[15]
    dflux_v18_v_17_dCsub = 1.0*RII_2*kf_RII_2__RII_C_2;  // dwdx[16]
    dflux_v20_v_19_dCsub = 1.0*RIIp_2*kf_RII_2__RII_C_2*xi_kf_RII_2__RII_C_2;  // dwdx[17]
    dflux_v40_v_39_dCsub = 1.0*H89*kf_H89;  // dwdx[18]
    dflux_v42_v_41_dCsub = 1.0*PDE*kf_PDE_Csub;  // dwdx[19]
    dflux_v24_v_23_dRIIp_C_2 = 1.0*cAMP*kf_cAMP;  // dwdx[25]
    dflux_v30_v_29_dRIIp_C_2 = 1.0*Rp8_Br_cAMPS*kf_cAMP*xi_b_Rp8_Br_cAMPS;  // dwdx[26]
    dflux_v32_v_31_dRIIp_C_2 = 1.0*Rp8_pCPT_cAMPS*kf_cAMP*xi_b_Rp8_pCPT_cAMPS;  // dwdx[27]
    dflux_v34_v_33_dRIIp_C_2 = 1.0*Rp_cAMPS*kf_cAMP*xi_b_Rp_cAMPS;  // dwdx[28]
    dflux_v36_v_35_dRIIp_C_2 = 1.0*Sp8_Br_cAMPS*kf_cAMP*xi_b_Sp8_Br_cAMPS;  // dwdx[29]
    dflux_v26_v_25_dRIIp_cAMP_C_2 = 1.0*KD_IBMX*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdx[31]
    dflux_v29_v_28_dRIIp_cAMP_2 = 1.0*KD_IBMX*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdx[34]
    dflux_v36_v_35_dSp8_Br_cAMPS = 1.0*RIIp_C_2*kf_cAMP*xi_b_Sp8_Br_cAMPS;  // dwdx[36]
    dflux_v30_v_29_dRp8_Br_cAMPS = 1.0*RIIp_C_2*kf_cAMP*xi_b_Rp8_Br_cAMPS;  // dwdx[38]
    dflux_v32_v_31_dRp8_pCPT_cAMPS = 1.0*RIIp_C_2*kf_cAMP*xi_b_Rp8_pCPT_cAMPS;  // dwdx[40]
    dflux_v34_v_33_dRp_cAMPS = 1.0*RIIp_C_2*kf_cAMP*xi_b_Rp_cAMPS;  // dwdx[42]
    dflux_v20_v_19_dRIIp_2 = 1.0*Csub*kf_RII_2__RII_C_2*xi_kf_RII_2__RII_C_2;  // dwdx[49]
    dflux_v18_v_17_dRII_2 = 1.0*Csub*kf_RII_2__RII_C_2;  // dwdx[51]
    dflux_v16_v_15_dIBMX = -1.0*KD_IBMX*cAMP*kdeg_cAMP_free*(PDE + pPDE*xi_pPDE)/std::pow(IBMX + KD_IBMX, 2);  // dwdx[53]
    dflux_v26_v_25_dIBMX = -1.0*KD_IBMX*RIIp_cAMP_C_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/std::pow(IBMX + KD_IBMX, 2);  // dwdx[54]
    dflux_v29_v_28_dIBMX = -1.0*KD_IBMX*RIIp_cAMP_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/std::pow(IBMX + KD_IBMX, 2);  // dwdx[55]
    dflux_v16_v_15_dPDE = 1.0*KD_IBMX*cAMP*kdeg_cAMP_free/(IBMX + KD_IBMX);  // dwdx[56]
    dflux_v26_v_25_dPDE = 1.0*KD_IBMX*RIIp_cAMP_C_2*kdeg_cAMP/(IBMX + KD_IBMX);  // dwdx[57]
    dflux_v29_v_28_dPDE = 1.0*KD_IBMX*RIIp_cAMP_2*kdeg_cAMP/(IBMX + KD_IBMX);  // dwdx[58]
    dflux_v42_v_41_dPDE = 1.0*Csub*kf_PDE_Csub;  // dwdx[59]
    dflux_v16_v_15_dpPDE = 1.0*KD_IBMX*cAMP*kdeg_cAMP_free*xi_pPDE/(IBMX + KD_IBMX);  // dwdx[60]
    dflux_v26_v_25_dpPDE = 1.0*KD_IBMX*RIIp_cAMP_C_2*kdeg_cAMP*xi_pPDE/(IBMX + KD_IBMX);  // dwdx[61]
    dflux_v29_v_28_dpPDE = 1.0*KD_IBMX*RIIp_cAMP_2*kdeg_cAMP*xi_pPDE/(IBMX + KD_IBMX);  // dwdx[62]
}

} // namespace model_Isensee_JCB2018
} // namespace amici
