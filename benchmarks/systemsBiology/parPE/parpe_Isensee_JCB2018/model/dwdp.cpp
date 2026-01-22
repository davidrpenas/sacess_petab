#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 70> dwdp_colptrs_Isensee_JCB2018_ = {
    0, 2, 3, 6, 13, 15, 16, 20, 22, 24, 26, 27, 28, 29, 31, 43, 44, 45, 46, 47, 48, 50, 51, 53, 55, 57, 60, 61, 62, 63, 64, 65, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67
};

void dwdp_colptrs_Isensee_JCB2018(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Isensee_JCB2018_));
}
} // namespace model_Isensee_JCB2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 67> dwdp_rowvals_Isensee_JCB2018_ = {
    14, 16, 48, 23, 33, 36, 32, 35, 38, 40, 42, 44, 46, 33, 36, 23, 13, 14, 15, 16, 47, 48, 25, 27, 24, 26, 29, 28, 30, 34, 45, 31, 32, 35, 37, 38, 39, 40, 41, 42, 43, 44, 46, 12, 9, 10, 11, 8, 21, 22, 22, 37, 38, 39, 40, 41, 42, 43, 44, 46, 27, 26, 38, 40, 42, 44, 46
};

void dwdp_rowvals_Isensee_JCB2018(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Isensee_JCB2018_));
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
#include "dwdp.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void dwdp_Isensee_JCB2018(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v7_v_6_dKD_Fsk = 1.0*AC_Fsk*kf_Fsk;  // dwdp[0]
    dflux_v9_v_8_dKD_Fsk = 1.0*kf_Fsk*pAC_Fsk;  // dwdp[1]
    dflux_v41_v_40_dKD_H89 = 1.0*Csub_H89*kf_H89;  // dwdp[2]
    dflux_v16_v_15_dKD_IBMX = -1.0*KD_IBMX*cAMP*kdeg_cAMP_free*(PDE + pPDE*xi_pPDE)/std::pow(IBMX + KD_IBMX, 2) + 1.0*cAMP*kdeg_cAMP_free*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdp[3]
    dflux_v26_v_25_dKD_IBMX = -1.0*KD_IBMX*RIIp_cAMP_C_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/std::pow(IBMX + KD_IBMX, 2) + 1.0*RIIp_cAMP_C_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdp[4]
    dflux_v29_v_28_dKD_IBMX = -1.0*KD_IBMX*RIIp_cAMP_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/std::pow(IBMX + KD_IBMX, 2) + 1.0*RIIp_cAMP_2*kdeg_cAMP*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdp[5]
    dflux_v25_v_24_dKD_cAMP = 1.0*RIIp_cAMP_C_2*kf_cAMP;  // dwdp[6]
    dflux_v28_v_27_dKD_cAMP = 1.0*RIIp_cAMP_2*kf_cAMP;  // dwdp[7]
    dflux_v31_v_30_dKD_cAMP = 1.0*RIIp_Rp8_Br_cAMPS_C_2*kf_cAMP*xi_KD_Rp8_Br_cAMPS*xi_b_Rp8_Br_cAMPS;  // dwdp[8]
    dflux_v33_v_32_dKD_cAMP = 1.0*RIIp_Rp8_pCPT_cAMPS_C_2*kf_cAMP*xi_KD_Rp8_pCPT_cAMPS*xi_b_Rp8_pCPT_cAMPS;  // dwdp[9]
    dflux_v35_v_34_dKD_cAMP = 1.0*RIIp_Rp_cAMPS_C_2*kf_cAMP*xi_KD_Rp_cAMPS*xi_b_Rp_cAMPS;  // dwdp[10]
    dflux_v37_v_36_dKD_cAMP = 1.0*RIIp_Sp8_Br_cAMPS_C_2*kf_cAMP*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // dwdp[11]
    dflux_v39_v_38_dKD_cAMP = 1.0*RIIp_Sp8_Br_cAMPS_2*kf_cAMP*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // dwdp[12]
    dflux_v26_v_25_dkdeg_cAMP = 1.0*KD_IBMX*RIIp_cAMP_C_2*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdp[13]
    dflux_v29_v_28_dkdeg_cAMP = 1.0*KD_IBMX*RIIp_cAMP_2*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdp[14]
    dflux_v16_v_15_dkdeg_cAMP_free = 1.0*KD_IBMX*cAMP*(PDE + pPDE*xi_pPDE)/(IBMX + KD_IBMX);  // dwdp[15]
    dflux_v6_v_5_dkf_Fsk = 1.0*AC*Fsk;  // dwdp[16]
    dflux_v7_v_6_dkf_Fsk = 1.0*AC_Fsk*KD_Fsk;  // dwdp[17]
    dflux_v8_v_7_dkf_Fsk = 1.0*Fsk*pAC;  // dwdp[18]
    dflux_v9_v_8_dkf_Fsk = 1.0*KD_Fsk*pAC_Fsk;  // dwdp[19]
    dflux_v40_v_39_dkf_H89 = 1.0*Csub*H89;  // dwdp[20]
    dflux_v41_v_40_dkf_H89 = 1.0*Csub_H89*KD_H89;  // dwdp[21]
    dflux_v18_v_17_dkf_RII_2__RII_C_2 = 1.0*Csub*RII_2;  // dwdp[22]
    dflux_v20_v_19_dkf_RII_2__RII_C_2 = 1.0*Csub*RIIp_2*xi_kf_RII_2__RII_C_2;  // dwdp[23]
    dflux_v17_v_16_dkf_RII_C_2__RII_2 = 1.0*RII_C_2;  // dwdp[24]
    dflux_v19_v_18_dkf_RII_C_2__RII_2 = 1.0*RIIp_C_2*xi_kf_RII_C_2__RII_2;  // dwdp[25]
    dflux_v22_v_21_dkf_RII_C_2__RIIp_C_2 = 1.0*RII_C_2;  // dwdp[26]
    dflux_v21_v_20_dkf_RIIp_2__RII_2 = 1.0*RIIp_2;  // dwdp[27]
    dflux_v23_v_22_dkf_RIIp_C_2__RII_C_2 = 1.0*RIIp_C_2;  // dwdp[28]
    dflux_v27_v_26_dkf_RIIp_cAMP_C_2__RIIp_2 = 1.0*RIIp_cAMP_C_2;  // dwdp[29]
    dflux_v38_v_37_dkf_RIIp_cAMP_C_2__RIIp_2 = 1.0*RIIp_Sp8_Br_cAMPS_C_2;  // dwdp[30]
    dflux_v24_v_23_dkf_cAMP = 1.0*RIIp_C_2*cAMP;  // dwdp[31]
    dflux_v25_v_24_dkf_cAMP = 1.0*KD_cAMP*RIIp_cAMP_C_2;  // dwdp[32]
    dflux_v28_v_27_dkf_cAMP = 1.0*KD_cAMP*RIIp_cAMP_2;  // dwdp[33]
    dflux_v30_v_29_dkf_cAMP = 1.0*RIIp_C_2*Rp8_Br_cAMPS*xi_b_Rp8_Br_cAMPS;  // dwdp[34]
    dflux_v31_v_30_dkf_cAMP = 1.0*KD_cAMP*RIIp_Rp8_Br_cAMPS_C_2*xi_KD_Rp8_Br_cAMPS*xi_b_Rp8_Br_cAMPS;  // dwdp[35]
    dflux_v32_v_31_dkf_cAMP = 1.0*RIIp_C_2*Rp8_pCPT_cAMPS*xi_b_Rp8_pCPT_cAMPS;  // dwdp[36]
    dflux_v33_v_32_dkf_cAMP = 1.0*KD_cAMP*RIIp_Rp8_pCPT_cAMPS_C_2*xi_KD_Rp8_pCPT_cAMPS*xi_b_Rp8_pCPT_cAMPS;  // dwdp[37]
    dflux_v34_v_33_dkf_cAMP = 1.0*RIIp_C_2*Rp_cAMPS*xi_b_Rp_cAMPS;  // dwdp[38]
    dflux_v35_v_34_dkf_cAMP = 1.0*KD_cAMP*RIIp_Rp_cAMPS_C_2*xi_KD_Rp_cAMPS*xi_b_Rp_cAMPS;  // dwdp[39]
    dflux_v36_v_35_dkf_cAMP = 1.0*RIIp_C_2*Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // dwdp[40]
    dflux_v37_v_36_dkf_cAMP = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_C_2*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // dwdp[41]
    dflux_v39_v_38_dkf_cAMP = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_2*xi_KD_Sp8_Br_cAMPS*xi_b_Sp8_Br_cAMPS;  // dwdp[42]
    dflux_v5_v_4_dki_IBMX = -1.0*IBMX + 1.0*IBMXex;  // dwdp[43]
    dflux_v2_v_1_dki_Rp8_Br_cAMPS_pAB = -1.0*Rp8_Br_cAMPS + 1.0*Rp8_Br_cAMPS_pAB*xi_i_Rp8_Br_cAMPS_pAB;  // dwdp[44]
    dflux_v3_v_2_dki_Rp8_pCPT_cAMPS_pAB = -1.0*Rp8_pCPT_cAMPS + 1.0*Rp8_pCPT_cAMPS_pAB*xi_i_Rp8_pCPT_cAMPS_pAB;  // dwdp[45]
    dflux_v4_v_3_dki_Rp_cAMPS_pAB = -1.0*Rp_cAMPS + 1.0*Rp_cAMPS_pAB*xi_i_Rp_cAMPS_pAB;  // dwdp[46]
    dflux_v1_v_0_dki_Sp8_Br_cAMPS_AM = -1.0*Sp8_Br_cAMPS + 1.0*Sp8_Br_cAMPS_AM*xi_i_Sp8_Br_cAMPS_AM;  // dwdp[47]
    dflux_v14_v_13_dks_AC_cAMP = 1.0*AC + 1.0*pAC*xi_pAC;  // dwdp[48]
    dflux_v15_v_14_dks_AC_cAMP = 1.0*xi_AC_cAMP_Fsk*(AC_Fsk + pAC_Fsk*xi_pAC);  // dwdp[49]
    dflux_v15_v_14_dxi_AC_cAMP_Fsk = 1.0*ks_AC_cAMP*(AC_Fsk + pAC_Fsk*xi_pAC);  // dwdp[50]
    dflux_v30_v_29_dxi_b_Rp8_Br_cAMPS = 1.0*RIIp_C_2*Rp8_Br_cAMPS*kf_cAMP;  // dwdp[51]
    dflux_v31_v_30_dxi_b_Rp8_Br_cAMPS = 1.0*KD_cAMP*RIIp_Rp8_Br_cAMPS_C_2*kf_cAMP*xi_KD_Rp8_Br_cAMPS;  // dwdp[52]
    dflux_v32_v_31_dxi_b_Rp8_pCPT_cAMPS = 1.0*RIIp_C_2*Rp8_pCPT_cAMPS*kf_cAMP;  // dwdp[53]
    dflux_v33_v_32_dxi_b_Rp8_pCPT_cAMPS = 1.0*KD_cAMP*RIIp_Rp8_pCPT_cAMPS_C_2*kf_cAMP*xi_KD_Rp8_pCPT_cAMPS;  // dwdp[54]
    dflux_v34_v_33_dxi_b_Rp_cAMPS = 1.0*RIIp_C_2*Rp_cAMPS*kf_cAMP;  // dwdp[55]
    dflux_v35_v_34_dxi_b_Rp_cAMPS = 1.0*KD_cAMP*RIIp_Rp_cAMPS_C_2*kf_cAMP*xi_KD_Rp_cAMPS;  // dwdp[56]
    dflux_v36_v_35_dxi_b_Sp8_Br_cAMPS = 1.0*RIIp_C_2*Sp8_Br_cAMPS*kf_cAMP;  // dwdp[57]
    dflux_v37_v_36_dxi_b_Sp8_Br_cAMPS = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_C_2*kf_cAMP*xi_KD_Sp8_Br_cAMPS;  // dwdp[58]
    dflux_v39_v_38_dxi_b_Sp8_Br_cAMPS = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_2*kf_cAMP*xi_KD_Sp8_Br_cAMPS;  // dwdp[59]
    dflux_v20_v_19_dxi_kf_RII_2__RII_C_2 = 1.0*Csub*RIIp_2*kf_RII_2__RII_C_2;  // dwdp[60]
    dflux_v19_v_18_dxi_kf_RII_C_2__RII_2 = 1.0*RIIp_C_2*kf_RII_C_2__RII_2;  // dwdp[61]
    dflux_v31_v_30_dxi_KD_Rp8_Br_cAMPS = 1.0*KD_cAMP*RIIp_Rp8_Br_cAMPS_C_2*kf_cAMP*xi_b_Rp8_Br_cAMPS;  // dwdp[62]
    dflux_v33_v_32_dxi_KD_Rp8_pCPT_cAMPS = 1.0*KD_cAMP*RIIp_Rp8_pCPT_cAMPS_C_2*kf_cAMP*xi_b_Rp8_pCPT_cAMPS;  // dwdp[63]
    dflux_v35_v_34_dxi_KD_Rp_cAMPS = 1.0*KD_cAMP*RIIp_Rp_cAMPS_C_2*kf_cAMP*xi_b_Rp_cAMPS;  // dwdp[64]
    dflux_v37_v_36_dxi_KD_Sp8_Br_cAMPS = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_C_2*kf_cAMP*xi_b_Sp8_Br_cAMPS;  // dwdp[65]
    dflux_v39_v_38_dxi_KD_Sp8_Br_cAMPS = 1.0*KD_cAMP*RIIp_Sp8_Br_cAMPS_2*kf_cAMP*xi_b_Sp8_Br_cAMPS;  // dwdp[66]
}

} // namespace model_Isensee_JCB2018
} // namespace amici
