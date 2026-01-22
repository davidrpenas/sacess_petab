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

void Jy_Raimundez_PCB2020(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EGFR_ID10_HS746T_FM_EGF_CET, 2)) + 0.5*std::pow(EGFR_ID10_HS746T_FM_EGF_CET - mEGFR_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_EGFR_ID10_HS746T_FM_EGF_CET, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EGFR_ID10_HS746T_FM_EGF, 2)) + 0.5*std::pow(EGFR_ID10_HS746T_FM_EGF - mEGFR_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_EGFR_ID10_HS746T_FM_EGF, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EGFR_ID12_MKN1_FM_HM, 2)) + 0.5*std::pow(EGFR_ID12_MKN1_FM_HM - mEGFR_ID12_MKN1_FM_HM, 2)/std::pow(sigma_EGFR_ID12_MKN1_FM_HM, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EGFR_ID3_HS746T_FM, 2)) + 0.5*std::pow(EGFR_ID3_HS746T_FM - mEGFR_ID3_HS746T_FM, 2)/std::pow(sigma_EGFR_ID3_HS746T_FM, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EGFR_ID3_MKN1_FM, 2)) + 0.5*std::pow(EGFR_ID3_MKN1_FM - mEGFR_ID3_MKN1_FM, 2)/std::pow(sigma_EGFR_ID3_MKN1_FM, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EGFR_ID4_HS746T_HM_3min, 2)) + 0.5*std::pow(EGFR_ID4_HS746T_HM_3min - mEGFR_ID4_HS746T_HM_3min, 2)/std::pow(sigma_EGFR_ID4_HS746T_HM_3min, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_EGFR_ID4_MKN1_HM_3min, 2)) + 0.5*std::pow(EGFR_ID4_MKN1_HM_3min - mEGFR_ID4_MKN1_HM_3min, 2)/std::pow(sigma_EGFR_ID4_MKN1_HM_3min, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MMET_ID10_HS746T_FM_EGF_CET, 2)) + 0.5*std::pow(MMET_ID10_HS746T_FM_EGF_CET - mMMET_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_MMET_ID10_HS746T_FM_EGF_CET, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MMET_ID10_HS746T_FM_EGF, 2)) + 0.5*std::pow(MMET_ID10_HS746T_FM_EGF - mMMET_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_MMET_ID10_HS746T_FM_EGF, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID12_MKN1_FM_HM, 2)) + 0.5*std::pow(-mpAKT_ID12_MKN1_FM_HM + pAKT_ID12_MKN1_FM_HM, 2)/std::pow(sigma_pAKT_ID12_MKN1_FM_HM, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID2_HS746T_HM, 2)) + 0.5*std::pow(-mpAKT_ID2_HS746T_HM + pAKT_ID2_HS746T_HM, 2)/std::pow(sigma_pAKT_ID2_HS746T_HM, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID2_MKN1_HM, 2)) + 0.5*std::pow(-mpAKT_ID2_MKN1_HM + pAKT_ID2_MKN1_HM, 2)/std::pow(sigma_pAKT_ID2_MKN1_HM, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID5_MKN1_FM_15min, 2)) + 0.5*std::pow(-mpAKT_ID5_MKN1_FM_15min + pAKT_ID5_MKN1_FM_15min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_15min, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID5_MKN1_FM_1min, 2)) + 0.5*std::pow(-mpAKT_ID5_MKN1_FM_1min + pAKT_ID5_MKN1_FM_1min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_1min, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID5_MKN1_FM_240min, 2)) + 0.5*std::pow(-mpAKT_ID5_MKN1_FM_240min + pAKT_ID5_MKN1_FM_240min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_240min, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID5_MKN1_FM_3min, 2)) + 0.5*std::pow(-mpAKT_ID5_MKN1_FM_3min + pAKT_ID5_MKN1_FM_3min, 2)/std::pow(sigma_pAKT_ID5_MKN1_FM_3min, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID7_HS746T_FM_15min, 2)) + 0.5*std::pow(-mpAKT_ID7_HS746T_FM_15min + pAKT_ID7_HS746T_FM_15min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_15min, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID7_HS746T_FM_1min, 2)) + 0.5*std::pow(-mpAKT_ID7_HS746T_FM_1min + pAKT_ID7_HS746T_FM_1min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_1min, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID7_HS746T_FM_240min, 2)) + 0.5*std::pow(-mpAKT_ID7_HS746T_FM_240min + pAKT_ID7_HS746T_FM_240min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_240min, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID7_HS746T_FM_3min, 2)) + 0.5*std::pow(-mpAKT_ID7_HS746T_FM_3min + pAKT_ID7_HS746T_FM_3min, 2)/std::pow(sigma_pAKT_ID7_HS746T_FM_3min, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID8_HS746T_HM_15min, 2)) + 0.5*std::pow(-mpAKT_ID8_HS746T_HM_15min + pAKT_ID8_HS746T_HM_15min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_15min, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID8_HS746T_HM_1min, 2)) + 0.5*std::pow(-mpAKT_ID8_HS746T_HM_1min + pAKT_ID8_HS746T_HM_1min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_1min, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID8_HS746T_HM_240min, 2)) + 0.5*std::pow(-mpAKT_ID8_HS746T_HM_240min + pAKT_ID8_HS746T_HM_240min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_240min, 2);
            break;
        case 23:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID8_HS746T_HM_3min, 2)) + 0.5*std::pow(-mpAKT_ID8_HS746T_HM_3min + pAKT_ID8_HS746T_HM_3min, 2)/std::pow(sigma_pAKT_ID8_HS746T_HM_3min, 2);
            break;
        case 24:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID9_HS746T_FM, 2)) + 0.5*std::pow(-mpAKT_ID9_HS746T_FM + pAKT_ID9_HS746T_FM, 2)/std::pow(sigma_pAKT_ID9_HS746T_FM, 2);
            break;
        case 25:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID9_HS746T_HM_30EGF, 2)) + 0.5*std::pow(-mpAKT_ID9_HS746T_HM_30EGF + pAKT_ID9_HS746T_HM_30EGF, 2)/std::pow(sigma_pAKT_ID9_HS746T_HM_30EGF, 2);
            break;
        case 26:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAKT_ID9_HS746T_HM_5EGF, 2)) + 0.5*std::pow(-mpAKT_ID9_HS746T_HM_5EGF + pAKT_ID9_HS746T_HM_5EGF, 2)/std::pow(sigma_pAKT_ID9_HS746T_HM_5EGF, 2);
            break;
        case 27:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID10_HS746T_FM_EGF_CET, 2)) + 0.5*std::pow(-mpEGFR_ID10_HS746T_FM_EGF_CET + pEGFR_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_pEGFR_ID10_HS746T_FM_EGF_CET, 2);
            break;
        case 28:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID10_HS746T_FM_EGF, 2)) + 0.5*std::pow(-mpEGFR_ID10_HS746T_FM_EGF + pEGFR_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_pEGFR_ID10_HS746T_FM_EGF, 2);
            break;
        case 29:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID12_MKN1_FM_HM, 2)) + 0.5*std::pow(-mpEGFR_ID12_MKN1_FM_HM + pEGFR_ID12_MKN1_FM_HM, 2)/std::pow(sigma_pEGFR_ID12_MKN1_FM_HM, 2);
            break;
        case 30:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID2_HS746T_HM, 2)) + 0.5*std::pow(-mpEGFR_ID2_HS746T_HM + pEGFR_ID2_HS746T_HM, 2)/std::pow(sigma_pEGFR_ID2_HS746T_HM, 2);
            break;
        case 31:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID2_MKN1_HM, 2)) + 0.5*std::pow(-mpEGFR_ID2_MKN1_HM + pEGFR_ID2_MKN1_HM, 2)/std::pow(sigma_pEGFR_ID2_MKN1_HM, 2);
            break;
        case 32:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID3_HS746T_FM, 2)) + 0.5*std::pow(-mpEGFR_ID3_HS746T_FM + pEGFR_ID3_HS746T_FM, 2)/std::pow(sigma_pEGFR_ID3_HS746T_FM, 2);
            break;
        case 33:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID3_MKN1_FM, 2)) + 0.5*std::pow(-mpEGFR_ID3_MKN1_FM + pEGFR_ID3_MKN1_FM, 2)/std::pow(sigma_pEGFR_ID3_MKN1_FM, 2);
            break;
        case 34:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID4_HS746T_HM_3min, 2)) + 0.5*std::pow(-mpEGFR_ID4_HS746T_HM_3min + pEGFR_ID4_HS746T_HM_3min, 2)/std::pow(sigma_pEGFR_ID4_HS746T_HM_3min, 2);
            break;
        case 35:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID4_MKN1_HM_3min, 2)) + 0.5*std::pow(-mpEGFR_ID4_MKN1_HM_3min + pEGFR_ID4_MKN1_HM_3min, 2)/std::pow(sigma_pEGFR_ID4_MKN1_HM_3min, 2);
            break;
        case 36:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID5_MKN1_FM_15min, 2)) + 0.5*std::pow(-mpEGFR_ID5_MKN1_FM_15min + pEGFR_ID5_MKN1_FM_15min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_15min, 2);
            break;
        case 37:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID5_MKN1_FM_1min, 2)) + 0.5*std::pow(-mpEGFR_ID5_MKN1_FM_1min + pEGFR_ID5_MKN1_FM_1min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_1min, 2);
            break;
        case 38:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID5_MKN1_FM_240min, 2)) + 0.5*std::pow(-mpEGFR_ID5_MKN1_FM_240min + pEGFR_ID5_MKN1_FM_240min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_240min, 2);
            break;
        case 39:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID5_MKN1_FM_3min, 2)) + 0.5*std::pow(-mpEGFR_ID5_MKN1_FM_3min + pEGFR_ID5_MKN1_FM_3min, 2)/std::pow(sigma_pEGFR_ID5_MKN1_FM_3min, 2);
            break;
        case 40:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID6_MKN1_HM_15min, 2)) + 0.5*std::pow(-mpEGFR_ID6_MKN1_HM_15min + pEGFR_ID6_MKN1_HM_15min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_15min, 2);
            break;
        case 41:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID6_MKN1_HM_1min, 2)) + 0.5*std::pow(-mpEGFR_ID6_MKN1_HM_1min + pEGFR_ID6_MKN1_HM_1min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_1min, 2);
            break;
        case 42:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID6_MKN1_HM_240min, 2)) + 0.5*std::pow(-mpEGFR_ID6_MKN1_HM_240min + pEGFR_ID6_MKN1_HM_240min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_240min, 2);
            break;
        case 43:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID6_MKN1_HM_3min, 2)) + 0.5*std::pow(-mpEGFR_ID6_MKN1_HM_3min + pEGFR_ID6_MKN1_HM_3min, 2)/std::pow(sigma_pEGFR_ID6_MKN1_HM_3min, 2);
            break;
        case 44:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID7_HS746T_FM_15min, 2)) + 0.5*std::pow(-mpEGFR_ID7_HS746T_FM_15min + pEGFR_ID7_HS746T_FM_15min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_15min, 2);
            break;
        case 45:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID7_HS746T_FM_1min, 2)) + 0.5*std::pow(-mpEGFR_ID7_HS746T_FM_1min + pEGFR_ID7_HS746T_FM_1min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_1min, 2);
            break;
        case 46:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID7_HS746T_FM_240min, 2)) + 0.5*std::pow(-mpEGFR_ID7_HS746T_FM_240min + pEGFR_ID7_HS746T_FM_240min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_240min, 2);
            break;
        case 47:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID7_HS746T_FM_3min, 2)) + 0.5*std::pow(-mpEGFR_ID7_HS746T_FM_3min + pEGFR_ID7_HS746T_FM_3min, 2)/std::pow(sigma_pEGFR_ID7_HS746T_FM_3min, 2);
            break;
        case 48:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID8_HS746T_HM_15min, 2)) + 0.5*std::pow(-mpEGFR_ID8_HS746T_HM_15min + pEGFR_ID8_HS746T_HM_15min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_15min, 2);
            break;
        case 49:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID8_HS746T_HM_1min, 2)) + 0.5*std::pow(-mpEGFR_ID8_HS746T_HM_1min + pEGFR_ID8_HS746T_HM_1min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_1min, 2);
            break;
        case 50:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID8_HS746T_HM_240min, 2)) + 0.5*std::pow(-mpEGFR_ID8_HS746T_HM_240min + pEGFR_ID8_HS746T_HM_240min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_240min, 2);
            break;
        case 51:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID8_HS746T_HM_3min, 2)) + 0.5*std::pow(-mpEGFR_ID8_HS746T_HM_3min + pEGFR_ID8_HS746T_HM_3min, 2)/std::pow(sigma_pEGFR_ID8_HS746T_HM_3min, 2);
            break;
        case 52:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID9_HS746T_FM, 2)) + 0.5*std::pow(-mpEGFR_ID9_HS746T_FM + pEGFR_ID9_HS746T_FM, 2)/std::pow(sigma_pEGFR_ID9_HS746T_FM, 2);
            break;
        case 53:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID9_HS746T_HM_30EGF, 2)) + 0.5*std::pow(-mpEGFR_ID9_HS746T_HM_30EGF + pEGFR_ID9_HS746T_HM_30EGF, 2)/std::pow(sigma_pEGFR_ID9_HS746T_HM_30EGF, 2);
            break;
        case 54:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID9_HS746T_HM_5EGF, 2)) + 0.5*std::pow(-mpEGFR_ID9_HS746T_HM_5EGF + pEGFR_ID9_HS746T_HM_5EGF, 2)/std::pow(sigma_pEGFR_ID9_HS746T_HM_5EGF, 2);
            break;
        case 55:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID9b_MKN1_HM_1EGF, 2)) + 0.5*std::pow(-mpEGFR_ID9b_MKN1_HM_1EGF + pEGFR_ID9b_MKN1_HM_1EGF, 2)/std::pow(sigma_pEGFR_ID9b_MKN1_HM_1EGF, 2);
            break;
        case 56:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_ID9b_MKN1_HM_5EGF, 2)) + 0.5*std::pow(-mpEGFR_ID9b_MKN1_HM_5EGF + pEGFR_ID9b_MKN1_HM_5EGF, 2)/std::pow(sigma_pEGFR_ID9b_MKN1_HM_5EGF, 2);
            break;
        case 57:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID10_HS746T_FM_EGF_CET, 2)) + 0.5*std::pow(-mpMAPK_ID10_HS746T_FM_EGF_CET + pMAPK_ID10_HS746T_FM_EGF_CET, 2)/std::pow(sigma_pMAPK_ID10_HS746T_FM_EGF_CET, 2);
            break;
        case 58:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID10_HS746T_FM_EGF, 2)) + 0.5*std::pow(-mpMAPK_ID10_HS746T_FM_EGF + pMAPK_ID10_HS746T_FM_EGF, 2)/std::pow(sigma_pMAPK_ID10_HS746T_FM_EGF, 2);
            break;
        case 59:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID12_MKN1_FM_HM, 2)) + 0.5*std::pow(-mpMAPK_ID12_MKN1_FM_HM + pMAPK_ID12_MKN1_FM_HM, 2)/std::pow(sigma_pMAPK_ID12_MKN1_FM_HM, 2);
            break;
        case 60:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID2_HS746T_HM, 2)) + 0.5*std::pow(-mpMAPK_ID2_HS746T_HM + pMAPK_ID2_HS746T_HM, 2)/std::pow(sigma_pMAPK_ID2_HS746T_HM, 2);
            break;
        case 61:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID2_MKN1_HM, 2)) + 0.5*std::pow(-mpMAPK_ID2_MKN1_HM + pMAPK_ID2_MKN1_HM, 2)/std::pow(sigma_pMAPK_ID2_MKN1_HM, 2);
            break;
        case 62:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID5_MKN1_FM_15min, 2)) + 0.5*std::pow(-mpMAPK_ID5_MKN1_FM_15min + pMAPK_ID5_MKN1_FM_15min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_15min, 2);
            break;
        case 63:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID5_MKN1_FM_1min, 2)) + 0.5*std::pow(-mpMAPK_ID5_MKN1_FM_1min + pMAPK_ID5_MKN1_FM_1min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_1min, 2);
            break;
        case 64:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID5_MKN1_FM_240min, 2)) + 0.5*std::pow(-mpMAPK_ID5_MKN1_FM_240min + pMAPK_ID5_MKN1_FM_240min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_240min, 2);
            break;
        case 65:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID5_MKN1_FM_3min, 2)) + 0.5*std::pow(-mpMAPK_ID5_MKN1_FM_3min + pMAPK_ID5_MKN1_FM_3min, 2)/std::pow(sigma_pMAPK_ID5_MKN1_FM_3min, 2);
            break;
        case 66:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID7_HS746T_FM_15min, 2)) + 0.5*std::pow(-mpMAPK_ID7_HS746T_FM_15min + pMAPK_ID7_HS746T_FM_15min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_15min, 2);
            break;
        case 67:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID7_HS746T_FM_1min, 2)) + 0.5*std::pow(-mpMAPK_ID7_HS746T_FM_1min + pMAPK_ID7_HS746T_FM_1min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_1min, 2);
            break;
        case 68:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID7_HS746T_FM_240min, 2)) + 0.5*std::pow(-mpMAPK_ID7_HS746T_FM_240min + pMAPK_ID7_HS746T_FM_240min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_240min, 2);
            break;
        case 69:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID7_HS746T_FM_3min, 2)) + 0.5*std::pow(-mpMAPK_ID7_HS746T_FM_3min + pMAPK_ID7_HS746T_FM_3min, 2)/std::pow(sigma_pMAPK_ID7_HS746T_FM_3min, 2);
            break;
        case 70:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID8_HS746T_HM_15min, 2)) + 0.5*std::pow(-mpMAPK_ID8_HS746T_HM_15min + pMAPK_ID8_HS746T_HM_15min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_15min, 2);
            break;
        case 71:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID8_HS746T_HM_1min, 2)) + 0.5*std::pow(-mpMAPK_ID8_HS746T_HM_1min + pMAPK_ID8_HS746T_HM_1min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_1min, 2);
            break;
        case 72:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID8_HS746T_HM_240min, 2)) + 0.5*std::pow(-mpMAPK_ID8_HS746T_HM_240min + pMAPK_ID8_HS746T_HM_240min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_240min, 2);
            break;
        case 73:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID8_HS746T_HM_3min, 2)) + 0.5*std::pow(-mpMAPK_ID8_HS746T_HM_3min + pMAPK_ID8_HS746T_HM_3min, 2)/std::pow(sigma_pMAPK_ID8_HS746T_HM_3min, 2);
            break;
        case 74:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID9_HS746T_FM, 2)) + 0.5*std::pow(-mpMAPK_ID9_HS746T_FM + pMAPK_ID9_HS746T_FM, 2)/std::pow(sigma_pMAPK_ID9_HS746T_FM, 2);
            break;
        case 75:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID9_HS746T_HM_30EGF, 2)) + 0.5*std::pow(-mpMAPK_ID9_HS746T_HM_30EGF + pMAPK_ID9_HS746T_HM_30EGF, 2)/std::pow(sigma_pMAPK_ID9_HS746T_HM_30EGF, 2);
            break;
        case 76:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID9_HS746T_HM_5EGF, 2)) + 0.5*std::pow(-mpMAPK_ID9_HS746T_HM_5EGF + pMAPK_ID9_HS746T_HM_5EGF, 2)/std::pow(sigma_pMAPK_ID9_HS746T_HM_5EGF, 2);
            break;
        case 77:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID9b_MKN1_HM_1EGF, 2)) + 0.5*std::pow(-mpMAPK_ID9b_MKN1_HM_1EGF + pMAPK_ID9b_MKN1_HM_1EGF, 2)/std::pow(sigma_pMAPK_ID9b_MKN1_HM_1EGF, 2);
            break;
        case 78:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pMAPK_ID9b_MKN1_HM_5EGF, 2)) + 0.5*std::pow(-mpMAPK_ID9b_MKN1_HM_5EGF + pMAPK_ID9b_MKN1_HM_5EGF, 2)/std::pow(sigma_pMAPK_ID9b_MKN1_HM_5EGF, 2);
            break;
        case 79:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_KD_EGFR_CET, 2)) + 0.5*std::pow(-mprior_KD_EGFR_CET + prior_KD_EGFR_CET, 2)/std::pow(sigma_prior_KD_EGFR_CET, 2);
            break;
        case 80:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_KD_EGFR_EGF, 2)) + 0.5*std::pow(-mprior_KD_EGFR_EGF + prior_KD_EGFR_EGF, 2)/std::pow(sigma_prior_KD_EGFR_EGF, 2);
            break;
        case 81:
            Jy[0] = std::log(2*sigma_prior_d_AKTtotal__MKN1_2_HS746T) + std::fabs(mprior_d_AKTtotal__MKN1_2_HS746T - prior_d_AKTtotal__MKN1_2_HS746T)/sigma_prior_d_AKTtotal__MKN1_2_HS746T;
            break;
        case 82:
            Jy[0] = std::log(2*sigma_prior_d_AKTtotal__fm_2_hm) + std::fabs(mprior_d_AKTtotal__fm_2_hm - prior_d_AKTtotal__fm_2_hm)/sigma_prior_d_AKTtotal__fm_2_hm;
            break;
        case 83:
            Jy[0] = std::log(2*sigma_prior_d_MAPKtotal__MKN1_2_HS746T) + std::fabs(mprior_d_MAPKtotal__MKN1_2_HS746T - prior_d_MAPKtotal__MKN1_2_HS746T)/sigma_prior_d_MAPKtotal__MKN1_2_HS746T;
            break;
        case 84:
            Jy[0] = std::log(2*sigma_prior_d_MAPKtotal__fm_2_hm) + std::fabs(mprior_d_MAPKtotal__fm_2_hm - prior_d_MAPKtotal__fm_2_hm)/sigma_prior_d_MAPKtotal__fm_2_hm;
            break;
        case 85:
            Jy[0] = std::log(2*sigma_prior_d_MPI3Ktotal__fm_2_hm) + std::fabs(mprior_d_MPI3Ktotal__fm_2_hm - prior_d_MPI3Ktotal__fm_2_hm)/sigma_prior_d_MPI3Ktotal__fm_2_hm;
            break;
        case 86:
            Jy[0] = std::log(2*sigma_prior_d_PI3Ktotal__fm_2_hm) + std::fabs(mprior_d_PI3Ktotal__fm_2_hm - prior_d_PI3Ktotal__fm_2_hm)/sigma_prior_d_PI3Ktotal__fm_2_hm;
            break;
        case 87:
            Jy[0] = std::log(2*sigma_prior_d_RAStotal__MKN1_2_HS746T) + std::fabs(mprior_d_RAStotal__MKN1_2_HS746T - prior_d_RAStotal__MKN1_2_HS746T)/sigma_prior_d_RAStotal__MKN1_2_HS746T;
            break;
        case 88:
            Jy[0] = std::log(2*sigma_prior_d_RAStotal__fm_2_hm) + std::fabs(mprior_d_RAStotal__fm_2_hm - prior_d_RAStotal__fm_2_hm)/sigma_prior_d_RAStotal__fm_2_hm;
            break;
        case 89:
            Jy[0] = std::log(2*sigma_prior_d_kdeg_membran__MKN1_2_HS746T) + std::fabs(mprior_d_kdeg_membran__MKN1_2_HS746T - prior_d_kdeg_membran__MKN1_2_HS746T)/sigma_prior_d_kdeg_membran__MKN1_2_HS746T;
            break;
        case 90:
            Jy[0] = std::log(2*sigma_prior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T) + std::fabs(mprior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T - prior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T)/sigma_prior_d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T;
            break;
        case 91:
            Jy[0] = std::log(2*sigma_prior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T) + std::fabs(mprior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T - prior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T)/sigma_prior_d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T;
            break;
        case 92:
            Jy[0] = std::log(2*sigma_prior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T) + std::fabs(mprior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T - prior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T)/sigma_prior_d_kimp_pEGFR_EGF_2__MKN1_2_HS746T;
            break;
        case 93:
            Jy[0] = std::log(2*sigma_prior_d_ksyn_EGFR__MKN1_2_HS746T) + std::fabs(mprior_d_ksyn_EGFR__MKN1_2_HS746T - prior_d_ksyn_EGFR__MKN1_2_HS746T)/sigma_prior_d_ksyn_EGFR__MKN1_2_HS746T;
            break;
        case 94:
            Jy[0] = std::log(2*sigma_prior_d_ksyn_EGFR__fm_2_hm) + std::fabs(mprior_d_ksyn_EGFR__fm_2_hm - prior_d_ksyn_EGFR__fm_2_hm)/sigma_prior_d_ksyn_EGFR__fm_2_hm;
            break;
        case 95:
            Jy[0] = std::log(2*sigma_prior_d_ksyn_MMET__fm_2_hm) + std::fabs(mprior_d_ksyn_MMET__fm_2_hm - prior_d_ksyn_MMET__fm_2_hm)/sigma_prior_d_ksyn_MMET__fm_2_hm;
            break;
        case 96:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_ki_AKT__MKN1, 2)) + 0.5*std::pow(-mprior_ki_AKT__MKN1 + prior_ki_AKT__MKN1, 2)/std::pow(sigma_prior_ki_AKT__MKN1, 2);
            break;
        case 97:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_kimp_pEGFR_EGF_2__MKN1, 2)) + 0.5*std::pow(-mprior_kimp_pEGFR_EGF_2__MKN1 + prior_kimp_pEGFR_EGF_2__MKN1, 2)/std::pow(sigma_prior_kimp_pEGFR_EGF_2__MKN1, 2);
            break;
        case 98:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_ka_PI3K_pMMET_2, 2)) + 0.5*std::pow(-mprior_xi_ka_PI3K_pMMET_2 + prior_xi_ka_PI3K_pMMET_2, 2)/std::pow(sigma_prior_xi_ka_PI3K_pMMET_2, 2);
            break;
        case 99:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_ka_PI3K_pMMET_pEGFR, 2)) + 0.5*std::pow(-mprior_xi_ka_PI3K_pMMET_pEGFR + prior_xi_ka_PI3K_pMMET_pEGFR, 2)/std::pow(sigma_prior_xi_ka_PI3K_pMMET_pEGFR, 2);
            break;
        case 100:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_ka_RAS_pMMET_2, 2)) + 0.5*std::pow(-mprior_xi_ka_RAS_pMMET_2 + prior_xi_ka_RAS_pMMET_2, 2)/std::pow(sigma_prior_xi_ka_RAS_pMMET_2, 2);
            break;
        case 101:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_ka_RAS_pMMET_pEGFR, 2)) + 0.5*std::pow(-mprior_xi_ka_RAS_pMMET_pEGFR + prior_xi_ka_RAS_pMMET_pEGFR, 2)/std::pow(sigma_prior_xi_ka_RAS_pMMET_pEGFR, 2);
            break;
        case 102:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kdeg_pMMET_2_i, 2)) + 0.5*std::pow(-mprior_xi_kdeg_pMMET_2_i + prior_xi_kdeg_pMMET_2_i, 2)/std::pow(sigma_prior_xi_kdeg_pMMET_2_i, 2);
            break;
        case 103:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kdeg_pMMET_pEGFR_i, 2)) + 0.5*std::pow(-mprior_xi_kdeg_pMMET_pEGFR_i + prior_xi_kdeg_pMMET_pEGFR_i, 2)/std::pow(sigma_prior_xi_kdeg_pMMET_pEGFR_i, 2);
            break;
        case 104:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kdim_MMET, 2)) + 0.5*std::pow(-mprior_xi_kdim_MMET + prior_xi_kdim_MMET, 2)/std::pow(sigma_prior_xi_kdim_MMET, 2);
            break;
        case 105:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kdim_MMET_EGFR, 2)) + 0.5*std::pow(-mprior_xi_kdim_MMET_EGFR + prior_xi_kdim_MMET_EGFR, 2)/std::pow(sigma_prior_xi_kdim_MMET_EGFR, 2);
            break;
        case 106:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kexp_pMMET_2_i, 2)) + 0.5*std::pow(-mprior_xi_kexp_pMMET_2_i + prior_xi_kexp_pMMET_2_i, 2)/std::pow(sigma_prior_xi_kexp_pMMET_2_i, 2);
            break;
        case 107:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kexp_pMMET_pEGFR_i, 2)) + 0.5*std::pow(-mprior_xi_kexp_pMMET_pEGFR_i + prior_xi_kexp_pMMET_pEGFR_i, 2)/std::pow(sigma_prior_xi_kexp_pMMET_pEGFR_i, 2);
            break;
        case 108:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_ki_MPI3K, 2)) + 0.5*std::pow(-mprior_xi_ki_MPI3K + prior_xi_ki_MPI3K, 2)/std::pow(sigma_prior_xi_ki_MPI3K, 2);
            break;
        case 109:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kimp_pMMET_2, 2)) + 0.5*std::pow(-mprior_xi_kimp_pMMET_2 + prior_xi_kimp_pMMET_2, 2)/std::pow(sigma_prior_xi_kimp_pMMET_2, 2);
            break;
        case 110:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kimp_pMMET_pEGFR, 2)) + 0.5*std::pow(-mprior_xi_kimp_pMMET_pEGFR + prior_xi_kimp_pMMET_pEGFR, 2)/std::pow(sigma_prior_xi_kimp_pMMET_pEGFR, 2);
            break;
    }
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
