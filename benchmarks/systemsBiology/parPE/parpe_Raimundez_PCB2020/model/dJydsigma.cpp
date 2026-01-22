#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void dJydsigma_Raimundez_PCB2020(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_EGFR_ID10_HS746T_FM_EGF_CET - 1.0*std::pow(EGFR_ID10_HS746T_FM_EGF_CET - mEGFR_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_EGFR_ID10_HS746T_FM_EGF_CET, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_EGFR_ID10_HS746T_FM_EGF - 1.0*std::pow(EGFR_ID10_HS746T_FM_EGF - mEGFR_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_EGFR_ID10_HS746T_FM_EGF, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_EGFR_ID12_MKN1_FM_HM - 1.0*std::pow(EGFR_ID12_MKN1_FM_HM - mEGFR_ID12_MKN1_FM_HM, 2)/std::pow(sigma_EGFR_ID12_MKN1_FM_HM, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_EGFR_ID3_HS746T_FM - 1.0*std::pow(EGFR_ID3_HS746T_FM - mEGFR_ID3_HS746T_FM, 2)/std::pow(sigma_EGFR_ID3_HS746T_FM, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_EGFR_ID3_MKN1_FM - 1.0*std::pow(EGFR_ID3_MKN1_FM - mEGFR_ID3_MKN1_FM, 2)/std::pow(sigma_EGFR_ID3_MKN1_FM, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_EGFR_ID4_HS746T_HM_3min - 1.0*std::pow(EGFR_ID4_HS746T_HM_3min - mEGFR_ID4_HS746T_HM_3min, 2)/std::pow(sigma_EGFR_ID4_HS746T_HM_3min, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_EGFR_ID4_MKN1_HM_3min - 1.0*std::pow(EGFR_ID4_MKN1_HM_3min - mEGFR_ID4_MKN1_HM_3min, 2)/std::pow(sigma_EGFR_ID4_MKN1_HM_3min, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_MMET_ID10_HS746T_FM_EGF_CET - 1.0*std::pow(MMET_ID10_HS746T_FM_EGF_CET - mMMET_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_MMET_ID10_HS746T_FM_EGF_CET, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_MMET_ID10_HS746T_FM_EGF - 1.0*std::pow(MMET_ID10_HS746T_FM_EGF - mMMET_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_MMET_ID10_HS746T_FM_EGF, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_pAKT_ID12_MKN1_FM_HM - 1.0*std::pow(-mpAKT_ID12_MKN1_FM_HM + pAKT_ID12_MKN1_FM_HM, 2)/std::pow(sigma_pAKT_ID12_MKN1_FM_HM, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_pAKT_ID2_HS746T_HM - 1.0*std::pow(-mpAKT_ID2_HS746T_HM + pAKT_ID2_HS746T_HM, 2)/std::pow(sigma_pAKT_ID2_HS746T_HM, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_pAKT_ID2_MKN1_HM - 1.0*std::pow(-mpAKT_ID2_MKN1_HM + pAKT_ID2_MKN1_HM, 2)/std::pow(sigma_pAKT_ID2_MKN1_HM, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_pAKT_ID5_MKN1_FM_15min - 1.0*std::pow(-mpAKT_ID5_MKN1_FM_15min + pAKT_ID5_MKN1_FM_15min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_15min, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_pAKT_ID5_MKN1_FM_1min - 1.0*std::pow(-mpAKT_ID5_MKN1_FM_1min + pAKT_ID5_MKN1_FM_1min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_1min, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_pAKT_ID5_MKN1_FM_240min - 1.0*std::pow(-mpAKT_ID5_MKN1_FM_240min + pAKT_ID5_MKN1_FM_240min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_240min, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_pAKT_ID5_MKN1_FM_3min - 1.0*std::pow(-mpAKT_ID5_MKN1_FM_3min + pAKT_ID5_MKN1_FM_3min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_3min, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_pAKT_ID7_HS746T_FM_15min - 1.0*std::pow(-mpAKT_ID7_HS746T_FM_15min + pAKT_ID7_HS746T_FM_15min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_15min, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_pAKT_ID7_HS746T_FM_1min - 1.0*std::pow(-mpAKT_ID7_HS746T_FM_1min + pAKT_ID7_HS746T_FM_1min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_1min, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_pAKT_ID7_HS746T_FM_240min - 1.0*std::pow(-mpAKT_ID7_HS746T_FM_240min + pAKT_ID7_HS746T_FM_240min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_240min, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_pAKT_ID7_HS746T_FM_3min - 1.0*std::pow(-mpAKT_ID7_HS746T_FM_3min + pAKT_ID7_HS746T_FM_3min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_3min, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_pAKT_ID8_HS746T_HM_15min - 1.0*std::pow(-mpAKT_ID8_HS746T_HM_15min + pAKT_ID8_HS746T_HM_15min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_15min, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_pAKT_ID8_HS746T_HM_1min - 1.0*std::pow(-mpAKT_ID8_HS746T_HM_1min + pAKT_ID8_HS746T_HM_1min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_1min, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_pAKT_ID8_HS746T_HM_240min - 1.0*std::pow(-mpAKT_ID8_HS746T_HM_240min + pAKT_ID8_HS746T_HM_240min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_240min, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_pAKT_ID8_HS746T_HM_3min - 1.0*std::pow(-mpAKT_ID8_HS746T_HM_3min + pAKT_ID8_HS746T_HM_3min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_3min, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_pAKT_ID9_HS746T_FM - 1.0*std::pow(-mpAKT_ID9_HS746T_FM + pAKT_ID9_HS746T_FM, 2)/std::pow(sigma_pAKT_ID9_HS746T_FM, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_pAKT_ID9_HS746T_HM_30EGF - 1.0*std::pow(-mpAKT_ID9_HS746T_HM_30EGF + pAKT_ID9_HS746T_HM_30EGF, 2)/std::pow(sigma_pAKT_ID9_HS746T_HM_30EGF, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_pAKT_ID9_HS746T_HM_5EGF - 1.0*std::pow(-mpAKT_ID9_HS746T_HM_5EGF + pAKT_ID9_HS746T_HM_5EGF, 2)/std::pow(sigma_pAKT_ID9_HS746T_HM_5EGF, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_pEGFR_ID10_HS746T_FM_EGF_CET - 1.0*std::pow(-mpEGFR_ID10_HS746T_FM_EGF_CET + pEGFR_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_pEGFR_ID10_HS746T_FM_EGF_CET, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_pEGFR_ID10_HS746T_FM_EGF - 1.0*std::pow(-mpEGFR_ID10_HS746T_FM_EGF + pEGFR_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_pEGFR_ID10_HS746T_FM_EGF, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_pEGFR_ID12_MKN1_FM_HM - 1.0*std::pow(-mpEGFR_ID12_MKN1_FM_HM + pEGFR_ID12_MKN1_FM_HM, 2)/std::pow(sigma_pEGFR_ID12_MKN1_FM_HM, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_pEGFR_ID2_HS746T_HM - 1.0*std::pow(-mpEGFR_ID2_HS746T_HM + pEGFR_ID2_HS746T_HM, 2)/std::pow(sigma_pEGFR_ID2_HS746T_HM, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_pEGFR_ID2_MKN1_HM - 1.0*std::pow(-mpEGFR_ID2_MKN1_HM + pEGFR_ID2_MKN1_HM, 2)/std::pow(sigma_pEGFR_ID2_MKN1_HM, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_pEGFR_ID3_HS746T_FM - 1.0*std::pow(-mpEGFR_ID3_HS746T_FM + pEGFR_ID3_HS746T_FM, 2)/std::pow(sigma_pEGFR_ID3_HS746T_FM, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_pEGFR_ID3_MKN1_FM - 1.0*std::pow(-mpEGFR_ID3_MKN1_FM + pEGFR_ID3_MKN1_FM, 2)/std::pow(sigma_pEGFR_ID3_MKN1_FM, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_pEGFR_ID4_HS746T_HM_3min - 1.0*std::pow(-mpEGFR_ID4_HS746T_HM_3min + pEGFR_ID4_HS746T_HM_3min, 2)/std::pow(sigma_pEGFR_ID4_HS746T_HM_3min, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_pEGFR_ID4_MKN1_HM_3min - 1.0*std::pow(-mpEGFR_ID4_MKN1_HM_3min + pEGFR_ID4_MKN1_HM_3min, 2)/std::pow(sigma_pEGFR_ID4_MKN1_HM_3min, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_pEGFR_ID5_MKN1_FM_15min - 1.0*std::pow(-mpEGFR_ID5_MKN1_FM_15min + pEGFR_ID5_MKN1_FM_15min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_15min, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_pEGFR_ID5_MKN1_FM_1min - 1.0*std::pow(-mpEGFR_ID5_MKN1_FM_1min + pEGFR_ID5_MKN1_FM_1min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_1min, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_pEGFR_ID5_MKN1_FM_240min - 1.0*std::pow(-mpEGFR_ID5_MKN1_FM_240min + pEGFR_ID5_MKN1_FM_240min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_240min, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_pEGFR_ID5_MKN1_FM_3min - 1.0*std::pow(-mpEGFR_ID5_MKN1_FM_3min + pEGFR_ID5_MKN1_FM_3min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_3min, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_pEGFR_ID6_MKN1_HM_15min - 1.0*std::pow(-mpEGFR_ID6_MKN1_HM_15min + pEGFR_ID6_MKN1_HM_15min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_15min, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_pEGFR_ID6_MKN1_HM_1min - 1.0*std::pow(-mpEGFR_ID6_MKN1_HM_1min + pEGFR_ID6_MKN1_HM_1min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_1min, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_pEGFR_ID6_MKN1_HM_240min - 1.0*std::pow(-mpEGFR_ID6_MKN1_HM_240min + pEGFR_ID6_MKN1_HM_240min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_240min, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_pEGFR_ID6_MKN1_HM_3min - 1.0*std::pow(-mpEGFR_ID6_MKN1_HM_3min + pEGFR_ID6_MKN1_HM_3min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_3min, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_pEGFR_ID7_HS746T_FM_15min - 1.0*std::pow(-mpEGFR_ID7_HS746T_FM_15min + pEGFR_ID7_HS746T_FM_15min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_15min, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_pEGFR_ID7_HS746T_FM_1min - 1.0*std::pow(-mpEGFR_ID7_HS746T_FM_1min + pEGFR_ID7_HS746T_FM_1min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_1min, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_pEGFR_ID7_HS746T_FM_240min - 1.0*std::pow(-mpEGFR_ID7_HS746T_FM_240min + pEGFR_ID7_HS746T_FM_240min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_240min, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_pEGFR_ID7_HS746T_FM_3min - 1.0*std::pow(-mpEGFR_ID7_HS746T_FM_3min + pEGFR_ID7_HS746T_FM_3min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_3min, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_pEGFR_ID8_HS746T_HM_15min - 1.0*std::pow(-mpEGFR_ID8_HS746T_HM_15min + pEGFR_ID8_HS746T_HM_15min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_15min, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_pEGFR_ID8_HS746T_HM_1min - 1.0*std::pow(-mpEGFR_ID8_HS746T_HM_1min + pEGFR_ID8_HS746T_HM_1min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_1min, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_pEGFR_ID8_HS746T_HM_240min - 1.0*std::pow(-mpEGFR_ID8_HS746T_HM_240min + pEGFR_ID8_HS746T_HM_240min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_240min, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_pEGFR_ID8_HS746T_HM_3min - 1.0*std::pow(-mpEGFR_ID8_HS746T_HM_3min + pEGFR_ID8_HS746T_HM_3min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_3min, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_pEGFR_ID9_HS746T_FM - 1.0*std::pow(-mpEGFR_ID9_HS746T_FM + pEGFR_ID9_HS746T_FM, 2)/std::pow(sigma_pEGFR_ID9_HS746T_FM, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_pEGFR_ID9_HS746T_HM_30EGF - 1.0*std::pow(-mpEGFR_ID9_HS746T_HM_30EGF + pEGFR_ID9_HS746T_HM_30EGF, 2)/std::pow(sigma_pEGFR_ID9_HS746T_HM_30EGF, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_pEGFR_ID9_HS746T_HM_5EGF - 1.0*std::pow(-mpEGFR_ID9_HS746T_HM_5EGF + pEGFR_ID9_HS746T_HM_5EGF, 2)/std::pow(sigma_pEGFR_ID9_HS746T_HM_5EGF, 3);
            break;
        case 55:
            dJydsigma[55] = 1.0/sigma_pEGFR_ID9b_MKN1_HM_1EGF - 1.0*std::pow(-mpEGFR_ID9b_MKN1_HM_1EGF + pEGFR_ID9b_MKN1_HM_1EGF, 2)/std::pow(sigma_pEGFR_ID9b_MKN1_HM_1EGF, 3);
            break;
        case 56:
            dJydsigma[56] = 1.0/sigma_pEGFR_ID9b_MKN1_HM_5EGF - 1.0*std::pow(-mpEGFR_ID9b_MKN1_HM_5EGF + pEGFR_ID9b_MKN1_HM_5EGF, 2)/std::pow(sigma_pEGFR_ID9b_MKN1_HM_5EGF, 3);
            break;
        case 57:
            dJydsigma[57] = 1.0/sigma_pMAPK_ID10_HS746T_FM_EGF_CET - 1.0*std::pow(-mpMAPK_ID10_HS746T_FM_EGF_CET + pMAPK_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_pMAPK_ID10_HS746T_FM_EGF_CET, 3);
            break;
        case 58:
            dJydsigma[58] = 1.0/sigma_pMAPK_ID10_HS746T_FM_EGF - 1.0*std::pow(-mpMAPK_ID10_HS746T_FM_EGF + pMAPK_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_pMAPK_ID10_HS746T_FM_EGF, 3);
            break;
        case 59:
            dJydsigma[59] = 1.0/sigma_pMAPK_ID12_MKN1_FM_HM - 1.0*std::pow(-mpMAPK_ID12_MKN1_FM_HM + pMAPK_ID12_MKN1_FM_HM, 2)/std::pow(sigma_pMAPK_ID12_MKN1_FM_HM, 3);
            break;
        case 60:
            dJydsigma[60] = 1.0/sigma_pMAPK_ID2_HS746T_HM - 1.0*std::pow(-mpMAPK_ID2_HS746T_HM + pMAPK_ID2_HS746T_HM, 2)/std::pow(sigma_pMAPK_ID2_HS746T_HM, 3);
            break;
        case 61:
            dJydsigma[61] = 1.0/sigma_pMAPK_ID2_MKN1_HM - 1.0*std::pow(-mpMAPK_ID2_MKN1_HM + pMAPK_ID2_MKN1_HM, 2)/std::pow(sigma_pMAPK_ID2_MKN1_HM, 3);
            break;
        case 62:
            dJydsigma[62] = 1.0/sigma_pMAPK_ID5_MKN1_FM_15min - 1.0*std::pow(-mpMAPK_ID5_MKN1_FM_15min + pMAPK_ID5_MKN1_FM_15min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_15min, 3);
            break;
        case 63:
            dJydsigma[63] = 1.0/sigma_pMAPK_ID5_MKN1_FM_1min - 1.0*std::pow(-mpMAPK_ID5_MKN1_FM_1min + pMAPK_ID5_MKN1_FM_1min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_1min, 3);
            break;
        case 64:
            dJydsigma[64] = 1.0/sigma_pMAPK_ID5_MKN1_FM_240min - 1.0*std::pow(-mpMAPK_ID5_MKN1_FM_240min + pMAPK_ID5_MKN1_FM_240min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_240min, 3);
            break;
        case 65:
            dJydsigma[65] = 1.0/sigma_pMAPK_ID5_MKN1_FM_3min - 1.0*std::pow(-mpMAPK_ID5_MKN1_FM_3min + pMAPK_ID5_MKN1_FM_3min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_3min, 3);
            break;
        case 66:
            dJydsigma[66] = 1.0/sigma_pMAPK_ID7_HS746T_FM_15min - 1.0*std::pow(-mpMAPK_ID7_HS746T_FM_15min + pMAPK_ID7_HS746T_FM_15min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_15min, 3);
            break;
        case 67:
            dJydsigma[67] = 1.0/sigma_pMAPK_ID7_HS746T_FM_1min - 1.0*std::pow(-mpMAPK_ID7_HS746T_FM_1min + pMAPK_ID7_HS746T_FM_1min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_1min, 3);
            break;
        case 68:
            dJydsigma[68] = 1.0/sigma_pMAPK_ID7_HS746T_FM_240min - 1.0*std::pow(-mpMAPK_ID7_HS746T_FM_240min + pMAPK_ID7_HS746T_FM_240min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_240min, 3);
            break;
        case 69:
            dJydsigma[69] = 1.0/sigma_pMAPK_ID7_HS746T_FM_3min - 1.0*std::pow(-mpMAPK_ID7_HS746T_FM_3min + pMAPK_ID7_HS746T_FM_3min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_3min, 3);
            break;
        case 70:
            dJydsigma[70] = 1.0/sigma_pMAPK_ID8_HS746T_HM_15min - 1.0*std::pow(-mpMAPK_ID8_HS746T_HM_15min + pMAPK_ID8_HS746T_HM_15min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_15min, 3);
            break;
        case 71:
            dJydsigma[71] = 1.0/sigma_pMAPK_ID8_HS746T_HM_1min - 1.0*std::pow(-mpMAPK_ID8_HS746T_HM_1min + pMAPK_ID8_HS746T_HM_1min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_1min, 3);
            break;
        case 72:
            dJydsigma[72] = 1.0/sigma_pMAPK_ID8_HS746T_HM_240min - 1.0*std::pow(-mpMAPK_ID8_HS746T_HM_240min + pMAPK_ID8_HS746T_HM_240min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_240min, 3);
            break;
        case 73:
            dJydsigma[73] = 1.0/sigma_pMAPK_ID8_HS746T_HM_3min - 1.0*std::pow(-mpMAPK_ID8_HS746T_HM_3min + pMAPK_ID8_HS746T_HM_3min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_3min, 3);
            break;
        case 74:
            dJydsigma[74] = 1.0/sigma_pMAPK_ID9_HS746T_FM - 1.0*std::pow(-mpMAPK_ID9_HS746T_FM + pMAPK_ID9_HS746T_FM, 2)/std::pow(sigma_pMAPK_ID9_HS746T_FM, 3);
            break;
        case 75:
            dJydsigma[75] = 1.0/sigma_pMAPK_ID9_HS746T_HM_30EGF - 1.0*std::pow(-mpMAPK_ID9_HS746T_HM_30EGF + pMAPK_ID9_HS746T_HM_30EGF, 2)/std::pow(sigma_pMAPK_ID9_HS746T_HM_30EGF, 3);
            break;
        case 76:
            dJydsigma[76] = 1.0/sigma_pMAPK_ID9_HS746T_HM_5EGF - 1.0*std::pow(-mpMAPK_ID9_HS746T_HM_5EGF + pMAPK_ID9_HS746T_HM_5EGF, 2)/std::pow(sigma_pMAPK_ID9_HS746T_HM_5EGF, 3);
            break;
        case 77:
            dJydsigma[77] = 1.0/sigma_pMAPK_ID9b_MKN1_HM_1EGF - 1.0*std::pow(-mpMAPK_ID9b_MKN1_HM_1EGF + pMAPK_ID9b_MKN1_HM_1EGF, 2)/std::pow(sigma_pMAPK_ID9b_MKN1_HM_1EGF, 3);
            break;
        case 78:
            dJydsigma[78] = 1.0/sigma_pMAPK_ID9b_MKN1_HM_5EGF - 1.0*std::pow(-mpMAPK_ID9b_MKN1_HM_5EGF + pMAPK_ID9b_MKN1_HM_5EGF, 2)/std::pow(sigma_pMAPK_ID9b_MKN1_HM_5EGF, 3);
            break;
        case 79:
            dJydsigma[79] = 1.0/sigma_prior_KD_EGFR_CET - 1.0*std::pow(-mprior_KD_EGFR_CET + prior_KD_EGFR_CET, 2)/std::pow(sigma_prior_KD_EGFR_CET, 3);
            break;
        case 80:
            dJydsigma[80] = 1.0/sigma_prior_KD_EGFR_EGF - 1.0*std::pow(-mprior_KD_EGFR_EGF + prior_KD_EGFR_EGF, 2)/std::pow(sigma_prior_KD_EGFR_EGF, 3);
            break;
        case 81:
            dJydsigma[81] = 1.0/sigma_prior_d_AKTtotal__MKN1_2_HS746T - std::fabs(mprior_d_AKTtotal__MKN1_2_HS746T - prior_d_AKTtotal__MKN1_2_HS746T)/std::pow(sigma_prior_d_AKTtotal__MKN1_2_HS746T, 2);
            break;
        case 82:
            dJydsigma[82] = 1.0/sigma_prior_d_AKTtotal__fm_2_hm - std::fabs(mprior_d_AKTtotal__fm_2_hm - prior_d_AKTtotal__fm_2_hm)/std::pow(sigma_prior_d_AKTtotal__fm_2_hm, 2);
            break;
        case 83:
            dJydsigma[83] = 1.0/sigma_prior_d_MAPKtotal__MKN1_2_HS746T - std::fabs(mprior_d_MAPKtotal__MKN1_2_HS746T - prior_d_MAPKtotal__MKN1_2_HS746T)/std::pow(sigma_prior_d_MAPKtotal__MKN1_2_HS746T, 2);
            break;
        case 84:
            dJydsigma[84] = 1.0/sigma_prior_d_MAPKtotal__fm_2_hm - std::fabs(mprior_d_MAPKtotal__fm_2_hm - prior_d_MAPKtotal__fm_2_hm)/std::pow(sigma_prior_d_MAPKtotal__fm_2_hm, 2);
            break;
        case 85:
            dJydsigma[85] = 1.0/sigma_prior_d_MPI3Ktotal__fm_2_hm - std::fabs(mprior_d_MPI3Ktotal__fm_2_hm - prior_d_MPI3Ktotal__fm_2_hm)/std::pow(sigma_prior_d_MPI3Ktotal__fm_2_hm, 2);
            break;
        case 86:
            dJydsigma[86] = 1.0/sigma_prior_d_PI3Ktotal__fm_2_hm - std::fabs(mprior_d_PI3Ktotal__fm_2_hm - prior_d_PI3Ktotal__fm_2_hm)/std::pow(sigma_prior_d_PI3Ktotal__fm_2_hm, 2);
            break;
        case 87:
            dJydsigma[87] = 1.0/sigma_prior_d_RAStotal__MKN1_2_HS746T - std::fabs(mprior_d_RAStotal__MKN1_2_HS746T - prior_d_RAStotal__MKN1_2_HS746T)/std::pow(sigma_prior_d_RAStotal__MKN1_2_HS746T, 2);
            break;
        case 88:
            dJydsigma[88] = 1.0/sigma_prior_d_RAStotal__fm_2_hm - std::fabs(mprior_d_RAStotal__fm_2_hm - prior_d_RAStotal__fm_2_hm)/std::pow(sigma_prior_d_RAStotal__fm_2_hm, 2);
            break;
        case 89:
            dJydsigma[89] = 1.0/sigma_prior_d_kdeg_membran__MKN1_2_HS746T - std::fabs(mprior_d_kdeg_membran__MKN1_2_HS746T - prior_d_kdeg_membran__MKN1_2_HS746T)/std::pow(sigma_prior_d_kdeg_membran__MKN1_2_HS746T, 2);
            break;
        case 90:
            dJydsigma[90] = 1.0/sigma_prior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T - std::fabs(mprior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T - prior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T)/std::pow(sigma_prior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T, 2);
            break;
        case 91:
            dJydsigma[91] = 1.0/sigma_prior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T - std::fabs(mprior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T - prior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T)/std::pow(sigma_prior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T, 2);
            break;
        case 92:
            dJydsigma[92] = 1.0/sigma_prior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T - std::fabs(mprior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T - prior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T)/std::pow(sigma_prior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T, 2);
            break;
        case 93:
            dJydsigma[93] = 1.0/sigma_prior_d_ksyn_EGFR__MKN1_2_HS746T - std::fabs(mprior_d_ksyn_EGFR__MKN1_2_HS746T - prior_d_ksyn_EGFR__MKN1_2_HS746T)/std::pow(sigma_prior_d_ksyn_EGFR__MKN1_2_HS746T, 2);
            break;
        case 94:
            dJydsigma[94] = 1.0/sigma_prior_d_ksyn_EGFR__fm_2_hm - std::fabs(mprior_d_ksyn_EGFR__fm_2_hm - prior_d_ksyn_EGFR__fm_2_hm)/std::pow(sigma_prior_d_ksyn_EGFR__fm_2_hm, 2);
            break;
        case 95:
            dJydsigma[95] = 1.0/sigma_prior_d_ksyn_MMET__fm_2_hm - std::fabs(mprior_d_ksyn_MMET__fm_2_hm - prior_d_ksyn_MMET__fm_2_hm)/std::pow(sigma_prior_d_ksyn_MMET__fm_2_hm, 2);
            break;
        case 96:
            dJydsigma[96] = 1.0/sigma_prior_ki_AKT__MKN1 - 1.0*std::pow(-mprior_ki_AKT__MKN1 + prior_ki_AKT__MKN1, 2)/std::pow(sigma_prior_ki_AKT__MKN1, 3);
            break;
        case 97:
            dJydsigma[97] = 1.0/sigma_prior_kimp_pEGFR_EGF_2__MKN1 - 1.0*std::pow(-mprior_kimp_pEGFR_EGF_2__MKN1 + prior_kimp_pEGFR_EGF_2__MKN1, 2)/std::pow(sigma_prior_kimp_pEGFR_EGF_2__MKN1, 3);
            break;
        case 98:
            dJydsigma[98] = 1.0/sigma_prior_xi_ka_PI3K_pMMET_2 - 1.0*std::pow(-mprior_xi_ka_PI3K_pMMET_2 + prior_xi_ka_PI3K_pMMET_2, 2)/std::pow(sigma_prior_xi_ka_PI3K_pMMET_2, 3);
            break;
        case 99:
            dJydsigma[99] = 1.0/sigma_prior_xi_ka_PI3K_pMMET_pEGFR - 1.0*std::pow(-mprior_xi_ka_PI3K_pMMET_pEGFR + prior_xi_ka_PI3K_pMMET_pEGFR, 2)/std::pow(sigma_prior_xi_ka_PI3K_pMMET_pEGFR, 3);
            break;
        case 100:
            dJydsigma[100] = 1.0/sigma_prior_xi_ka_RAS_pMMET_2 - 1.0*std::pow(-mprior_xi_ka_RAS_pMMET_2 + prior_xi_ka_RAS_pMMET_2, 2)/std::pow(sigma_prior_xi_ka_RAS_pMMET_2, 3);
            break;
        case 101:
            dJydsigma[101] = 1.0/sigma_prior_xi_ka_RAS_pMMET_pEGFR - 1.0*std::pow(-mprior_xi_ka_RAS_pMMET_pEGFR + prior_xi_ka_RAS_pMMET_pEGFR, 2)/std::pow(sigma_prior_xi_ka_RAS_pMMET_pEGFR, 3);
            break;
        case 102:
            dJydsigma[102] = 1.0/sigma_prior_xi_kdeg_pMMET_2_i - 1.0*std::pow(-mprior_xi_kdeg_pMMET_2_i + prior_xi_kdeg_pMMET_2_i, 2)/std::pow(sigma_prior_xi_kdeg_pMMET_2_i, 3);
            break;
        case 103:
            dJydsigma[103] = 1.0/sigma_prior_xi_kdeg_pMMET_pEGFR_i - 1.0*std::pow(-mprior_xi_kdeg_pMMET_pEGFR_i + prior_xi_kdeg_pMMET_pEGFR_i, 2)/std::pow(sigma_prior_xi_kdeg_pMMET_pEGFR_i, 3);
            break;
        case 104:
            dJydsigma[104] = 1.0/sigma_prior_xi_kdim_MMET - 1.0*std::pow(-mprior_xi_kdim_MMET + prior_xi_kdim_MMET, 2)/std::pow(sigma_prior_xi_kdim_MMET, 3);
            break;
        case 105:
            dJydsigma[105] = 1.0/sigma_prior_xi_kdim_MMET_EGFR - 1.0*std::pow(-mprior_xi_kdim_MMET_EGFR + prior_xi_kdim_MMET_EGFR, 2)/std::pow(sigma_prior_xi_kdim_MMET_EGFR, 3);
            break;
        case 106:
            dJydsigma[106] = 1.0/sigma_prior_xi_kexp_pMMET_2_i - 1.0*std::pow(-mprior_xi_kexp_pMMET_2_i + prior_xi_kexp_pMMET_2_i, 2)/std::pow(sigma_prior_xi_kexp_pMMET_2_i, 3);
            break;
        case 107:
            dJydsigma[107] = 1.0/sigma_prior_xi_kexp_pMMET_pEGFR_i - 1.0*std::pow(-mprior_xi_kexp_pMMET_pEGFR_i + prior_xi_kexp_pMMET_pEGFR_i, 2)/std::pow(sigma_prior_xi_kexp_pMMET_pEGFR_i, 3);
            break;
        case 108:
            dJydsigma[108] = 1.0/sigma_prior_xi_ki_MPI3K - 1.0*std::pow(-mprior_xi_ki_MPI3K + prior_xi_ki_MPI3K, 2)/std::pow(sigma_prior_xi_ki_MPI3K, 3);
            break;
        case 109:
            dJydsigma[109] = 1.0/sigma_prior_xi_kimp_pMMET_2 - 1.0*std::pow(-mprior_xi_kimp_pMMET_2 + prior_xi_kimp_pMMET_2, 2)/std::pow(sigma_prior_xi_kimp_pMMET_2, 3);
            break;
        case 110:
            dJydsigma[110] = 1.0/sigma_prior_xi_kimp_pMMET_pEGFR - 1.0*std::pow(-mprior_xi_kimp_pMMET_pEGFR + prior_xi_kimp_pMMET_pEGFR, 2)/std::pow(sigma_prior_xi_kimp_pMMET_pEGFR, 3);
            break;
    }
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
