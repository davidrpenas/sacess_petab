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
namespace model_Raimundez_PCB2020 {

void y_Raimundez_PCB2020(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = observableParameter1_EGFR_ID10_HS746T_FM_EGF_CET*(EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[1] = observableParameter1_EGFR_ID10_HS746T_FM_EGF*(EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[2] = observableParameter1_EGFR_ID12_MKN1_FM_HM*(EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[3] = observableParameter1_EGFR_ID3_HS746T_FM*(EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[4] = observableParameter1_EGFR_ID3_MKN1_FM*(EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[5] = observableParameter1_EGFR_ID4_HS746T_HM_3min*(EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[6] = observableParameter1_EGFR_ID4_MKN1_HM_3min*(EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[7] = observableParameter1_MMET_ID10_HS746T_FM_EGF_CET*(EGFR_MMET_METinh + MMET + 2*MMET_2 + MMET_EGFR + MMET_METinh + 2*MMET_METinh_2 + 2*MMET_MMET_METinh + 2*pMMET_2 + 2*pMMET_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[8] = observableParameter1_MMET_ID10_HS746T_FM_EGF*(EGFR_MMET_METinh + MMET + 2*MMET_2 + MMET_EGFR + MMET_METinh + 2*MMET_METinh_2 + 2*MMET_MMET_METinh + 2*pMMET_2 + 2*pMMET_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[9] = observableParameter1_pAKT_ID12_MKN1_FM_HM*pAKT;
    y[10] = observableParameter1_pAKT_ID2_HS746T_HM*pAKT;
    y[11] = observableParameter1_pAKT_ID2_MKN1_HM*pAKT;
    y[12] = observableParameter1_pAKT_ID5_MKN1_FM_15min*pAKT;
    y[13] = observableParameter1_pAKT_ID5_MKN1_FM_1min*pAKT;
    y[14] = observableParameter1_pAKT_ID5_MKN1_FM_240min*pAKT;
    y[15] = observableParameter1_pAKT_ID5_MKN1_FM_3min*pAKT;
    y[16] = observableParameter1_pAKT_ID7_HS746T_FM_15min*pAKT;
    y[17] = observableParameter1_pAKT_ID7_HS746T_FM_1min*pAKT;
    y[18] = observableParameter1_pAKT_ID7_HS746T_FM_240min*pAKT;
    y[19] = observableParameter1_pAKT_ID7_HS746T_FM_3min*pAKT;
    y[20] = observableParameter1_pAKT_ID8_HS746T_HM_15min*pAKT;
    y[21] = observableParameter1_pAKT_ID8_HS746T_HM_1min*pAKT;
    y[22] = observableParameter1_pAKT_ID8_HS746T_HM_240min*pAKT;
    y[23] = observableParameter1_pAKT_ID8_HS746T_HM_3min*pAKT;
    y[24] = observableParameter1_pAKT_ID9_HS746T_FM*pAKT;
    y[25] = observableParameter1_pAKT_ID9_HS746T_HM_30EGF*pAKT;
    y[26] = observableParameter1_pAKT_ID9_HS746T_HM_5EGF*pAKT;
    y[27] = observableParameter1_pEGFR_ID10_HS746T_FM_EGF_CET*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[28] = observableParameter1_pEGFR_ID10_HS746T_FM_EGF*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[29] = observableParameter1_pEGFR_ID12_MKN1_FM_HM*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[30] = observableParameter1_pEGFR_ID2_HS746T_HM*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[31] = observableParameter1_pEGFR_ID2_MKN1_HM*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[32] = observableParameter1_pEGFR_ID3_HS746T_FM*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[33] = observableParameter1_pEGFR_ID3_MKN1_FM*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[34] = observableParameter1_pEGFR_ID4_HS746T_HM_3min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[35] = observableParameter1_pEGFR_ID4_MKN1_HM_3min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[36] = observableParameter1_pEGFR_ID5_MKN1_FM_15min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[37] = observableParameter1_pEGFR_ID5_MKN1_FM_1min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[38] = observableParameter1_pEGFR_ID5_MKN1_FM_240min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[39] = observableParameter1_pEGFR_ID5_MKN1_FM_3min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[40] = observableParameter1_pEGFR_ID6_MKN1_HM_15min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[41] = observableParameter1_pEGFR_ID6_MKN1_HM_1min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[42] = observableParameter1_pEGFR_ID6_MKN1_HM_240min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[43] = observableParameter1_pEGFR_ID6_MKN1_HM_3min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[44] = observableParameter1_pEGFR_ID7_HS746T_FM_15min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[45] = observableParameter1_pEGFR_ID7_HS746T_FM_1min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[46] = observableParameter1_pEGFR_ID7_HS746T_FM_240min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[47] = observableParameter1_pEGFR_ID7_HS746T_FM_3min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[48] = observableParameter1_pEGFR_ID8_HS746T_HM_15min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[49] = observableParameter1_pEGFR_ID8_HS746T_HM_1min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[50] = observableParameter1_pEGFR_ID8_HS746T_HM_240min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[51] = observableParameter1_pEGFR_ID8_HS746T_HM_3min*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[52] = observableParameter1_pEGFR_ID9_HS746T_FM*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[53] = observableParameter1_pEGFR_ID9_HS746T_HM_30EGF*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[54] = observableParameter1_pEGFR_ID9_HS746T_HM_5EGF*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[55] = observableParameter1_pEGFR_ID9b_MKN1_HM_1EGF*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[56] = observableParameter1_pEGFR_ID9b_MKN1_HM_5EGF*(2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i);
    y[57] = observableParameter1_pMAPK_ID10_HS746T_FM_EGF_CET*pMAPK;
    y[58] = observableParameter1_pMAPK_ID10_HS746T_FM_EGF*pMAPK;
    y[59] = observableParameter1_pMAPK_ID12_MKN1_FM_HM*pMAPK;
    y[60] = observableParameter1_pMAPK_ID2_HS746T_HM*pMAPK;
    y[61] = observableParameter1_pMAPK_ID2_MKN1_HM*pMAPK;
    y[62] = observableParameter1_pMAPK_ID5_MKN1_FM_15min*pMAPK;
    y[63] = observableParameter1_pMAPK_ID5_MKN1_FM_1min*pMAPK;
    y[64] = observableParameter1_pMAPK_ID5_MKN1_FM_240min*pMAPK;
    y[65] = observableParameter1_pMAPK_ID5_MKN1_FM_3min*pMAPK;
    y[66] = observableParameter1_pMAPK_ID7_HS746T_FM_15min*pMAPK;
    y[67] = observableParameter1_pMAPK_ID7_HS746T_FM_1min*pMAPK;
    y[68] = observableParameter1_pMAPK_ID7_HS746T_FM_240min*pMAPK;
    y[69] = observableParameter1_pMAPK_ID7_HS746T_FM_3min*pMAPK;
    y[70] = observableParameter1_pMAPK_ID8_HS746T_HM_15min*pMAPK;
    y[71] = observableParameter1_pMAPK_ID8_HS746T_HM_1min*pMAPK;
    y[72] = observableParameter1_pMAPK_ID8_HS746T_HM_240min*pMAPK;
    y[73] = observableParameter1_pMAPK_ID8_HS746T_HM_3min*pMAPK;
    y[74] = observableParameter1_pMAPK_ID9_HS746T_FM*pMAPK;
    y[75] = observableParameter1_pMAPK_ID9_HS746T_HM_30EGF*pMAPK;
    y[76] = observableParameter1_pMAPK_ID9_HS746T_HM_5EGF*pMAPK;
    y[77] = observableParameter1_pMAPK_ID9b_MKN1_HM_1EGF*pMAPK;
    y[78] = observableParameter1_pMAPK_ID9b_MKN1_HM_5EGF*pMAPK;
    y[79] = std::log(KD_EGFR_CET)/M_LN10;
    y[80] = std::log(KD_EGFR_EGF)/M_LN10;
    y[81] = d_AKTtotal__MKN1_2_HS746T;
    y[82] = d_AKTtotal__fm_2_hm;
    y[83] = d_MAPKtotal__MKN1_2_HS746T;
    y[84] = d_MAPKtotal__fm_2_hm;
    y[85] = d_MPI3Ktotal__fm_2_hm;
    y[86] = d_PI3Ktotal__fm_2_hm;
    y[87] = d_RAStotal__MKN1_2_HS746T;
    y[88] = d_RAStotal__fm_2_hm;
    y[89] = d_kdeg_membran__MKN1_2_HS746T;
    y[90] = d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T;
    y[91] = d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T;
    y[92] = d_kimp_pEGFR_EGF_2__MKN1_2_HS746T;
    y[93] = d_ksyn_EGFR__MKN1_2_HS746T;
    y[94] = d_ksyn_EGFR__fm_2_hm;
    y[95] = d_ksyn_MMET__fm_2_hm;
    y[96] = std::log(ki_AKT__MKN1)/M_LN10;
    y[97] = std::log(kimp_pEGFR_EGF_2__MKN1)/M_LN10;
    y[98] = std::log(xi_ka_PI3K_pMMET_2)/M_LN10;
    y[99] = std::log(xi_ka_PI3K_pMMET_pEGFR)/M_LN10;
    y[100] = std::log(xi_ka_RAS_pMMET_2)/M_LN10;
    y[101] = std::log(xi_ka_RAS_pMMET_pEGFR)/M_LN10;
    y[102] = std::log(xi_kdeg_pMMET_2_i)/M_LN10;
    y[103] = std::log(xi_kdeg_pMMET_pEGFR_i)/M_LN10;
    y[104] = std::log(xi_kdim_MMET)/M_LN10;
    y[105] = std::log(xi_kdim_MMET_EGFR)/M_LN10;
    y[106] = std::log(xi_kexp_pMMET_2_i)/M_LN10;
    y[107] = std::log(xi_kexp_pMMET_pEGFR_i)/M_LN10;
    y[108] = std::log(xi_ki_MPI3K)/M_LN10;
    y[109] = std::log(xi_kimp_pMMET_2)/M_LN10;
    y[110] = std::log(xi_kimp_pMMET_pEGFR)/M_LN10;
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
