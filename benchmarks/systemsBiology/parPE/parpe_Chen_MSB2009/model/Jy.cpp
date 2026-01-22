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

void Jy_Chen_MSB2009(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__0_1__model1_data1 - mAKT_PP__AKT_t__0_1__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data1, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__0_1__model1_data2 - mAKT_PP__AKT_t__0_1__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data2, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__0_1__model1_data3 - mAKT_PP__AKT_t__0_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data3, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data4, 2)) + 0.5*std::pow(AKT_PP__AKT_t__0_1__model1_data4 - mAKT_PP__AKT_t__0_1__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__0_1__model1_data4, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__0_36__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__0_36__model1_data1 - mAKT_PP__AKT_t__0_36__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__0_36__model1_data1, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__0_69__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__0_69__model1_data1 - mAKT_PP__AKT_t__0_69__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__0_69__model1_data1, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__0_98__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__0_98__model1_data3 - mAKT_PP__AKT_t__0_98__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__0_98__model1_data3, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__1_5__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__1_5__model1_data1 - mAKT_PP__AKT_t__1_5__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__1_5__model1_data1, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__1_6__model1_data4, 2)) + 0.5*std::pow(AKT_PP__AKT_t__1_6__model1_data4 - mAKT_PP__AKT_t__1_6__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__1_6__model1_data4, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__1_8__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__1_8__model1_data2 - mAKT_PP__AKT_t__1_8__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__1_8__model1_data2, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__2_1__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__2_1__model1_data1 - mAKT_PP__AKT_t__2_1__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__2_1__model1_data1, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__2_3__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__2_3__model1_data1 - mAKT_PP__AKT_t__2_3__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__2_3__model1_data1, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__2_3__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__2_3__model1_data2 - mAKT_PP__AKT_t__2_3__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_3__model1_data2, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__2_6__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__2_6__model1_data2 - mAKT_PP__AKT_t__2_6__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_6__model1_data2, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__2_8__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__2_8__model1_data2 - mAKT_PP__AKT_t__2_8__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_8__model1_data2, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__2_9__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__2_9__model1_data2 - mAKT_PP__AKT_t__2_9__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__2_9__model1_data2, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__3_1__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__3_1__model1_data3 - mAKT_PP__AKT_t__3_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__3_1__model1_data3, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__3_2__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__3_2__model1_data2 - mAKT_PP__AKT_t__3_2__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__3_2__model1_data2, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__3_5__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__3_5__model1_data1 - mAKT_PP__AKT_t__3_5__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__3_5__model1_data1, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__3_6__model1_data4, 2)) + 0.5*std::pow(AKT_PP__AKT_t__3_6__model1_data4 - mAKT_PP__AKT_t__3_6__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__3_6__model1_data4, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__4_5__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__4_5__model1_data3 - mAKT_PP__AKT_t__4_5__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__4_5__model1_data3, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__4_6__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__4_6__model1_data2 - mAKT_PP__AKT_t__4_6__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__4_6__model1_data2, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__4_6__model1_data4, 2)) + 0.5*std::pow(AKT_PP__AKT_t__4_6__model1_data4 - mAKT_PP__AKT_t__4_6__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__4_6__model1_data4, 2);
            break;
        case 23:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__4_7__model1_data2, 2)) + 0.5*std::pow(AKT_PP__AKT_t__4_7__model1_data2 - mAKT_PP__AKT_t__4_7__model1_data2, 2)/std::pow(sigma_AKT_PP__AKT_t__4_7__model1_data2, 2);
            break;
        case 24:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__5_1__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__5_1__model1_data3 - mAKT_PP__AKT_t__5_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__5_1__model1_data3, 2);
            break;
        case 25:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__5_6__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__5_6__model1_data3 - mAKT_PP__AKT_t__5_6__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__5_6__model1_data3, 2);
            break;
        case 26:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__7_4__model1_data1, 2)) + 0.5*std::pow(AKT_PP__AKT_t__7_4__model1_data1 - mAKT_PP__AKT_t__7_4__model1_data1, 2)/std::pow(sigma_AKT_PP__AKT_t__7_4__model1_data1, 2);
            break;
        case 27:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__7_7__model1_data4, 2)) + 0.5*std::pow(AKT_PP__AKT_t__7_7__model1_data4 - mAKT_PP__AKT_t__7_7__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__7_7__model1_data4, 2);
            break;
        case 28:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__9_7__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__9_7__model1_data3 - mAKT_PP__AKT_t__9_7__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__9_7__model1_data3, 2);
            break;
        case 29:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__11_1__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__11_1__model1_data3 - mAKT_PP__AKT_t__11_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__11_1__model1_data3, 2);
            break;
        case 30:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__12_1__model1_data4, 2)) + 0.5*std::pow(AKT_PP__AKT_t__12_1__model1_data4 - mAKT_PP__AKT_t__12_1__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__12_1__model1_data4, 2);
            break;
        case 31:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__13_1__model1_data4, 2)) + 0.5*std::pow(AKT_PP__AKT_t__13_1__model1_data4 - mAKT_PP__AKT_t__13_1__model1_data4, 2)/std::pow(sigma_AKT_PP__AKT_t__13_1__model1_data4, 2);
            break;
        case 32:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_AKT_PP__AKT_t__14_1__model1_data3, 2)) + 0.5*std::pow(AKT_PP__AKT_t__14_1__model1_data3 - mAKT_PP__AKT_t__14_1__model1_data3, 2)/std::pow(sigma_AKT_PP__AKT_t__14_1__model1_data3, 2);
            break;
        case 33:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data1 - mERB_B1_P_tot__EGFR_t__0_1__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data1, 2);
            break;
        case 34:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data2, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data2 - mERB_B1_P_tot__EGFR_t__0_1__model1_data2, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data2, 2);
            break;
        case 35:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data3, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data3 - mERB_B1_P_tot__EGFR_t__0_1__model1_data3, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data3, 2);
            break;
        case 36:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data4, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__0_1__model1_data4 - mERB_B1_P_tot__EGFR_t__0_1__model1_data4, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data4, 2);
            break;
        case 37:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__2_4__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__2_4__model1_data1 - mERB_B1_P_tot__EGFR_t__2_4__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__2_4__model1_data1, 2);
            break;
        case 38:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__2_5__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__2_5__model1_data1 - mERB_B1_P_tot__EGFR_t__2_5__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__2_5__model1_data1, 2);
            break;
        case 39:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__4_7__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__4_7__model1_data1 - mERB_B1_P_tot__EGFR_t__4_7__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__4_7__model1_data1, 2);
            break;
        case 40:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__6_8__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__6_8__model1_data1 - mERB_B1_P_tot__EGFR_t__6_8__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__6_8__model1_data1, 2);
            break;
        case 41:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__7_4__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__7_4__model1_data1 - mERB_B1_P_tot__EGFR_t__7_4__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__7_4__model1_data1, 2);
            break;
        case 42:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__8_9__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__8_9__model1_data1 - mERB_B1_P_tot__EGFR_t__8_9__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__8_9__model1_data1, 2);
            break;
        case 43:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERB_B1_P_tot__EGFR_t__13_1__model1_data1, 2)) + 0.5*std::pow(ERB_B1_P_tot__EGFR_t__13_1__model1_data1 - mERB_B1_P_tot__EGFR_t__13_1__model1_data1, 2)/std::pow(sigma_ERB_B1_P_tot__EGFR_t__13_1__model1_data1, 2);
            break;
        case 44:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__0_1__model1_data1 - mERK_PP__ERK_t__0_1__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data1, 2);
            break;
        case 45:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__0_1__model1_data2 - mERK_PP__ERK_t__0_1__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data2, 2);
            break;
        case 46:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__0_1__model1_data3 - mERK_PP__ERK_t__0_1__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data3, 2);
            break;
        case 47:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__0_1__model1_data4 - mERK_PP__ERK_t__0_1__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__0_1__model1_data4, 2);
            break;
        case 48:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__0_46__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__0_46__model1_data1 - mERK_PP__ERK_t__0_46__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__0_46__model1_data1, 2);
            break;
        case 49:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__0_94__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__0_94__model1_data2 - mERK_PP__ERK_t__0_94__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__0_94__model1_data2, 2);
            break;
        case 50:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__0_94__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__0_94__model1_data3 - mERK_PP__ERK_t__0_94__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__0_94__model1_data3, 2);
            break;
        case 51:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_04__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_04__model1_data4 - mERK_PP__ERK_t__1_04__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_04__model1_data4, 2);
            break;
        case 52:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_1__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_1__model1_data4 - mERK_PP__ERK_t__1_1__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_1__model1_data4, 2);
            break;
        case 53:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_2__model1_data1 - mERK_PP__ERK_t__1_2__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data1, 2);
            break;
        case 54:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_2__model1_data2 - mERK_PP__ERK_t__1_2__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data2, 2);
            break;
        case 55:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_2__model1_data3 - mERK_PP__ERK_t__1_2__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__1_2__model1_data3, 2);
            break;
        case 56:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_5__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_5__model1_data1 - mERK_PP__ERK_t__1_5__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__1_5__model1_data1, 2);
            break;
        case 57:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_5__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_5__model1_data4 - mERK_PP__ERK_t__1_5__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_5__model1_data4, 2);
            break;
        case 58:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_7__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_7__model1_data4 - mERK_PP__ERK_t__1_7__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__1_7__model1_data4, 2);
            break;
        case 59:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_9__model1_data1 - mERK_PP__ERK_t__1_9__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data1, 2);
            break;
        case 60:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_9__model1_data2 - mERK_PP__ERK_t__1_9__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data2, 2);
            break;
        case 61:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__1_9__model1_data3 - mERK_PP__ERK_t__1_9__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__1_9__model1_data3, 2);
            break;
        case 62:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__2_2__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__2_2__model1_data4 - mERK_PP__ERK_t__2_2__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__2_2__model1_data4, 2);
            break;
        case 63:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__2_3__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__2_3__model1_data2 - mERK_PP__ERK_t__2_3__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__2_3__model1_data2, 2);
            break;
        case 64:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__2_3__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__2_3__model1_data3 - mERK_PP__ERK_t__2_3__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__2_3__model1_data3, 2);
            break;
        case 65:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__2_7__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__2_7__model1_data2 - mERK_PP__ERK_t__2_7__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__2_7__model1_data2, 2);
            break;
        case 66:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__2_7__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__2_7__model1_data3 - mERK_PP__ERK_t__2_7__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__2_7__model1_data3, 2);
            break;
        case 67:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__2_8__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__2_8__model1_data1 - mERK_PP__ERK_t__2_8__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__2_8__model1_data1, 2);
            break;
        case 68:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__2_9__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__2_9__model1_data1 - mERK_PP__ERK_t__2_9__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__2_9__model1_data1, 2);
            break;
        case 69:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__3_0__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__3_0__model1_data4 - mERK_PP__ERK_t__3_0__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__3_0__model1_data4, 2);
            break;
        case 70:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__3_1__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__3_1__model1_data1 - mERK_PP__ERK_t__3_1__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__3_1__model1_data1, 2);
            break;
        case 71:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__3_4__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__3_4__model1_data2 - mERK_PP__ERK_t__3_4__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__3_4__model1_data2, 2);
            break;
        case 72:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__3_4__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__3_4__model1_data3 - mERK_PP__ERK_t__3_4__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__3_4__model1_data3, 2);
            break;
        case 73:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__3_7__model1_data1, 2)) + 0.5*std::pow(ERK_PP__ERK_t__3_7__model1_data1 - mERK_PP__ERK_t__3_7__model1_data1, 2)/std::pow(sigma_ERK_PP__ERK_t__3_7__model1_data1, 2);
            break;
        case 74:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__4_0__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__4_0__model1_data2 - mERK_PP__ERK_t__4_0__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__4_0__model1_data2, 2);
            break;
        case 75:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__4_0__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__4_0__model1_data3 - mERK_PP__ERK_t__4_0__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__4_0__model1_data3, 2);
            break;
        case 76:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__4_3__model1_data2 - mERK_PP__ERK_t__4_3__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data2, 2);
            break;
        case 77:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__4_3__model1_data3 - mERK_PP__ERK_t__4_3__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data3, 2);
            break;
        case 78:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__4_3__model1_data4 - mERK_PP__ERK_t__4_3__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__4_3__model1_data4, 2);
            break;
        case 79:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__4_7__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__4_7__model1_data4 - mERK_PP__ERK_t__4_7__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__4_7__model1_data4, 2);
            break;
        case 80:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__6_8__model1_data2, 2)) + 0.5*std::pow(ERK_PP__ERK_t__6_8__model1_data2 - mERK_PP__ERK_t__6_8__model1_data2, 2)/std::pow(sigma_ERK_PP__ERK_t__6_8__model1_data2, 2);
            break;
        case 81:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__6_8__model1_data3, 2)) + 0.5*std::pow(ERK_PP__ERK_t__6_8__model1_data3 - mERK_PP__ERK_t__6_8__model1_data3, 2)/std::pow(sigma_ERK_PP__ERK_t__6_8__model1_data3, 2);
            break;
        case 82:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ERK_PP__ERK_t__8_2__model1_data4, 2)) + 0.5*std::pow(ERK_PP__ERK_t__8_2__model1_data4 - mERK_PP__ERK_t__8_2__model1_data4, 2)/std::pow(sigma_ERK_PP__ERK_t__8_2__model1_data4, 2);
            break;
    }
}

} // namespace model_Chen_MSB2009
} // namespace amici
