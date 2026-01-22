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
namespace model_Chen_MSB2009 {

void dJydsigma_Chen_MSB2009(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_AKT_PP__AKT_t__0_1__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__0_1__model1_data1 - mAKT_PP__AKT_t__0_1__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_AKT_PP__AKT_t__0_1__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__0_1__model1_data2 - mAKT_PP__AKT_t__0_1__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_AKT_PP__AKT_t__0_1__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__0_1__model1_data3 - mAKT_PP__AKT_t__0_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data3, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_AKT_PP__AKT_t__0_1__model1_data4 - 1.0*std::pow(AKT_PP__AKT_t__0_1__model1_data4 - mAKT_PP__AKT_t__0_1__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data4, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_AKT_PP__AKT_t__0_36__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__0_36__model1_data1 - mAKT_PP__AKT_t__0_36__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__0_36__model1_data1, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_AKT_PP__AKT_t__0_69__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__0_69__model1_data1 - mAKT_PP__AKT_t__0_69__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__0_69__model1_data1, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_AKT_PP__AKT_t__0_98__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__0_98__model1_data3 - mAKT_PP__AKT_t__0_98__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__0_98__model1_data3, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_AKT_PP__AKT_t__1_5__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__1_5__model1_data1 - mAKT_PP__AKT_t__1_5__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__1_5__model1_data1, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_AKT_PP__AKT_t__1_6__model1_data4 - 1.0*std::pow(AKT_PP__AKT_t__1_6__model1_data4 - mAKT_PP__AKT_t__1_6__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__1_6__model1_data4, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_AKT_PP__AKT_t__1_8__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__1_8__model1_data2 - mAKT_PP__AKT_t__1_8__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__1_8__model1_data2, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_AKT_PP__AKT_t__2_1__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__2_1__model1_data1 - mAKT_PP__AKT_t__2_1__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__2_1__model1_data1, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_AKT_PP__AKT_t__2_3__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__2_3__model1_data1 - mAKT_PP__AKT_t__2_3__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__2_3__model1_data1, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_AKT_PP__AKT_t__2_3__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__2_3__model1_data2 - mAKT_PP__AKT_t__2_3__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_3__model1_data2, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_AKT_PP__AKT_t__2_6__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__2_6__model1_data2 - mAKT_PP__AKT_t__2_6__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_6__model1_data2, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_AKT_PP__AKT_t__2_8__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__2_8__model1_data2 - mAKT_PP__AKT_t__2_8__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_8__model1_data2, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_AKT_PP__AKT_t__2_9__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__2_9__model1_data2 - mAKT_PP__AKT_t__2_9__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_9__model1_data2, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_AKT_PP__AKT_t__3_1__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__3_1__model1_data3 - mAKT_PP__AKT_t__3_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__3_1__model1_data3, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_AKT_PP__AKT_t__3_2__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__3_2__model1_data2 - mAKT_PP__AKT_t__3_2__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__3_2__model1_data2, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_AKT_PP__AKT_t__3_5__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__3_5__model1_data1 - mAKT_PP__AKT_t__3_5__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__3_5__model1_data1, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_AKT_PP__AKT_t__3_6__model1_data4 - 1.0*std::pow(AKT_PP__AKT_t__3_6__model1_data4 - mAKT_PP__AKT_t__3_6__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__3_6__model1_data4, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_AKT_PP__AKT_t__4_5__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__4_5__model1_data3 - mAKT_PP__AKT_t__4_5__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__4_5__model1_data3, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_AKT_PP__AKT_t__4_6__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__4_6__model1_data2 - mAKT_PP__AKT_t__4_6__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__4_6__model1_data2, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_AKT_PP__AKT_t__4_6__model1_data4 - 1.0*std::pow(AKT_PP__AKT_t__4_6__model1_data4 - mAKT_PP__AKT_t__4_6__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__4_6__model1_data4, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_AKT_PP__AKT_t__4_7__model1_data2 - 1.0*std::pow(AKT_PP__AKT_t__4_7__model1_data2 - mAKT_PP__AKT_t__4_7__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__4_7__model1_data2, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_AKT_PP__AKT_t__5_1__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__5_1__model1_data3 - mAKT_PP__AKT_t__5_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__5_1__model1_data3, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_AKT_PP__AKT_t__5_6__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__5_6__model1_data3 - mAKT_PP__AKT_t__5_6__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__5_6__model1_data3, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_AKT_PP__AKT_t__7_4__model1_data1 - 1.0*std::pow(AKT_PP__AKT_t__7_4__model1_data1 - mAKT_PP__AKT_t__7_4__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__7_4__model1_data1, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_AKT_PP__AKT_t__7_7__model1_data4 - 1.0*std::pow(AKT_PP__AKT_t__7_7__model1_data4 - mAKT_PP__AKT_t__7_7__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__7_7__model1_data4, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_AKT_PP__AKT_t__9_7__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__9_7__model1_data3 - mAKT_PP__AKT_t__9_7__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__9_7__model1_data3, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_AKT_PP__AKT_t__11_1__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__11_1__model1_data3 - mAKT_PP__AKT_t__11_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__11_1__model1_data3, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_AKT_PP__AKT_t__12_1__model1_data4 - 1.0*std::pow(AKT_PP__AKT_t__12_1__model1_data4 - mAKT_PP__AKT_t__12_1__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__12_1__model1_data4, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_AKT_PP__AKT_t__13_1__model1_data4 - 1.0*std::pow(AKT_PP__AKT_t__13_1__model1_data4 - mAKT_PP__AKT_t__13_1__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__13_1__model1_data4, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_AKT_PP__AKT_t__14_1__model1_data3 - 1.0*std::pow(AKT_PP__AKT_t__14_1__model1_data3 - mAKT_PP__AKT_t__14_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__14_1__model1_data3, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data1 - mERB_B1_P_tot__EGFR_t__0_1__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data1, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data2 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data2 - mERB_B1_P_tot__EGFR_t__0_1__model1_data2, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data2, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data3 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data3 - mERB_B1_P_tot__EGFR_t__0_1__model1_data3, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data3, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data4 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data4 - mERB_B1_P_tot__EGFR_t__0_1__model1_data4, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data4, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__2_4__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__2_4__model1_data1 - mERB_B1_P_tot__EGFR_t__2_4__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__2_4__model1_data1, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__2_5__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__2_5__model1_data1 - mERB_B1_P_tot__EGFR_t__2_5__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__2_5__model1_data1, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__4_7__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__4_7__model1_data1 - mERB_B1_P_tot__EGFR_t__4_7__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__4_7__model1_data1, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__6_8__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__6_8__model1_data1 - mERB_B1_P_tot__EGFR_t__6_8__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__6_8__model1_data1, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__7_4__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__7_4__model1_data1 - mERB_B1_P_tot__EGFR_t__7_4__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__7_4__model1_data1, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__8_9__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__8_9__model1_data1 - mERB_B1_P_tot__EGFR_t__8_9__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__8_9__model1_data1, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_ERB_B1_P_tot__EGFR_t__13_1__model1_data1 - 1.0*std::pow(ERB_B1_P_tot__EGFR_t__13_1__model1_data1 - mERB_B1_P_tot__EGFR_t__13_1__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__13_1__model1_data1, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_ERK_PP__ERK_t__0_1__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__0_1__model1_data1 - mERK_PP__ERK_t__0_1__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data1, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_ERK_PP__ERK_t__0_1__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__0_1__model1_data2 - mERK_PP__ERK_t__0_1__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data2, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_ERK_PP__ERK_t__0_1__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__0_1__model1_data3 - mERK_PP__ERK_t__0_1__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data3, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_ERK_PP__ERK_t__0_1__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__0_1__model1_data4 - mERK_PP__ERK_t__0_1__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data4, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_ERK_PP__ERK_t__0_46__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__0_46__model1_data1 - mERK_PP__ERK_t__0_46__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__0_46__model1_data1, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_ERK_PP__ERK_t__0_94__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__0_94__model1_data2 - mERK_PP__ERK_t__0_94__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__0_94__model1_data2, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_ERK_PP__ERK_t__0_94__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__0_94__model1_data3 - mERK_PP__ERK_t__0_94__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__0_94__model1_data3, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_ERK_PP__ERK_t__1_04__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__1_04__model1_data4 - mERK_PP__ERK_t__1_04__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_04__model1_data4, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_ERK_PP__ERK_t__1_1__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__1_1__model1_data4 - mERK_PP__ERK_t__1_1__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_1__model1_data4, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_ERK_PP__ERK_t__1_2__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__1_2__model1_data1 - mERK_PP__ERK_t__1_2__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data1, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_ERK_PP__ERK_t__1_2__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__1_2__model1_data2 - mERK_PP__ERK_t__1_2__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data2, 3);
            break;
        case 55:
            dJydsigma[55] = 1.0/sigma_ERK_PP__ERK_t__1_2__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__1_2__model1_data3 - mERK_PP__ERK_t__1_2__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data3, 3);
            break;
        case 56:
            dJydsigma[56] = 1.0/sigma_ERK_PP__ERK_t__1_5__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__1_5__model1_data1 - mERK_PP__ERK_t__1_5__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__1_5__model1_data1, 3);
            break;
        case 57:
            dJydsigma[57] = 1.0/sigma_ERK_PP__ERK_t__1_5__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__1_5__model1_data4 - mERK_PP__ERK_t__1_5__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_5__model1_data4, 3);
            break;
        case 58:
            dJydsigma[58] = 1.0/sigma_ERK_PP__ERK_t__1_7__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__1_7__model1_data4 - mERK_PP__ERK_t__1_7__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_7__model1_data4, 3);
            break;
        case 59:
            dJydsigma[59] = 1.0/sigma_ERK_PP__ERK_t__1_9__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__1_9__model1_data1 - mERK_PP__ERK_t__1_9__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data1, 3);
            break;
        case 60:
            dJydsigma[60] = 1.0/sigma_ERK_PP__ERK_t__1_9__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__1_9__model1_data2 - mERK_PP__ERK_t__1_9__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data2, 3);
            break;
        case 61:
            dJydsigma[61] = 1.0/sigma_ERK_PP__ERK_t__1_9__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__1_9__model1_data3 - mERK_PP__ERK_t__1_9__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data3, 3);
            break;
        case 62:
            dJydsigma[62] = 1.0/sigma_ERK_PP__ERK_t__2_2__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__2_2__model1_data4 - mERK_PP__ERK_t__2_2__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__2_2__model1_data4, 3);
            break;
        case 63:
            dJydsigma[63] = 1.0/sigma_ERK_PP__ERK_t__2_3__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__2_3__model1_data2 - mERK_PP__ERK_t__2_3__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__2_3__model1_data2, 3);
            break;
        case 64:
            dJydsigma[64] = 1.0/sigma_ERK_PP__ERK_t__2_3__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__2_3__model1_data3 - mERK_PP__ERK_t__2_3__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__2_3__model1_data3, 3);
            break;
        case 65:
            dJydsigma[65] = 1.0/sigma_ERK_PP__ERK_t__2_7__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__2_7__model1_data2 - mERK_PP__ERK_t__2_7__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__2_7__model1_data2, 3);
            break;
        case 66:
            dJydsigma[66] = 1.0/sigma_ERK_PP__ERK_t__2_7__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__2_7__model1_data3 - mERK_PP__ERK_t__2_7__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__2_7__model1_data3, 3);
            break;
        case 67:
            dJydsigma[67] = 1.0/sigma_ERK_PP__ERK_t__2_8__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__2_8__model1_data1 - mERK_PP__ERK_t__2_8__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__2_8__model1_data1, 3);
            break;
        case 68:
            dJydsigma[68] = 1.0/sigma_ERK_PP__ERK_t__2_9__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__2_9__model1_data1 - mERK_PP__ERK_t__2_9__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__2_9__model1_data1, 3);
            break;
        case 69:
            dJydsigma[69] = 1.0/sigma_ERK_PP__ERK_t__3_0__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__3_0__model1_data4 - mERK_PP__ERK_t__3_0__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__3_0__model1_data4, 3);
            break;
        case 70:
            dJydsigma[70] = 1.0/sigma_ERK_PP__ERK_t__3_1__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__3_1__model1_data1 - mERK_PP__ERK_t__3_1__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__3_1__model1_data1, 3);
            break;
        case 71:
            dJydsigma[71] = 1.0/sigma_ERK_PP__ERK_t__3_4__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__3_4__model1_data2 - mERK_PP__ERK_t__3_4__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__3_4__model1_data2, 3);
            break;
        case 72:
            dJydsigma[72] = 1.0/sigma_ERK_PP__ERK_t__3_4__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__3_4__model1_data3 - mERK_PP__ERK_t__3_4__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__3_4__model1_data3, 3);
            break;
        case 73:
            dJydsigma[73] = 1.0/sigma_ERK_PP__ERK_t__3_7__model1_data1 - 1.0*std::pow(ERK_PP__ERK_t__3_7__model1_data1 - mERK_PP__ERK_t__3_7__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__3_7__model1_data1, 3);
            break;
        case 74:
            dJydsigma[74] = 1.0/sigma_ERK_PP__ERK_t__4_0__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__4_0__model1_data2 - mERK_PP__ERK_t__4_0__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__4_0__model1_data2, 3);
            break;
        case 75:
            dJydsigma[75] = 1.0/sigma_ERK_PP__ERK_t__4_0__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__4_0__model1_data3 - mERK_PP__ERK_t__4_0__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__4_0__model1_data3, 3);
            break;
        case 76:
            dJydsigma[76] = 1.0/sigma_ERK_PP__ERK_t__4_3__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__4_3__model1_data2 - mERK_PP__ERK_t__4_3__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data2, 3);
            break;
        case 77:
            dJydsigma[77] = 1.0/sigma_ERK_PP__ERK_t__4_3__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__4_3__model1_data3 - mERK_PP__ERK_t__4_3__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data3, 3);
            break;
        case 78:
            dJydsigma[78] = 1.0/sigma_ERK_PP__ERK_t__4_3__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__4_3__model1_data4 - mERK_PP__ERK_t__4_3__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data4, 3);
            break;
        case 79:
            dJydsigma[79] = 1.0/sigma_ERK_PP__ERK_t__4_7__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__4_7__model1_data4 - mERK_PP__ERK_t__4_7__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__4_7__model1_data4, 3);
            break;
        case 80:
            dJydsigma[80] = 1.0/sigma_ERK_PP__ERK_t__6_8__model1_data2 - 1.0*std::pow(ERK_PP__ERK_t__6_8__model1_data2 - mERK_PP__ERK_t__6_8__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__6_8__model1_data2, 3);
            break;
        case 81:
            dJydsigma[81] = 1.0/sigma_ERK_PP__ERK_t__6_8__model1_data3 - 1.0*std::pow(ERK_PP__ERK_t__6_8__model1_data3 - mERK_PP__ERK_t__6_8__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__6_8__model1_data3, 3);
            break;
        case 82:
            dJydsigma[82] = 1.0/sigma_ERK_PP__ERK_t__8_2__model1_data4 - 1.0*std::pow(ERK_PP__ERK_t__8_2__model1_data4 - mERK_PP__ERK_t__8_2__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__8_2__model1_data4, 3);
            break;
    }
}

} // namespace model_Chen_MSB2009
} // namespace amici
