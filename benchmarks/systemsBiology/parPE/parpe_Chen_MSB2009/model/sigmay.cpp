#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Chen_MSB2009 {

void sigmay_Chen_MSB2009(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_AKT_PP__AKT_t__0_1__model1_data1 = noiseParameter1_AKT_PP__AKT_t__0_1__model1_data1;  // sigmay[0]
    sigma_AKT_PP__AKT_t__0_1__model1_data2 = noiseParameter1_AKT_PP__AKT_t__0_1__model1_data2;  // sigmay[1]
    sigma_AKT_PP__AKT_t__0_1__model1_data3 = noiseParameter1_AKT_PP__AKT_t__0_1__model1_data3;  // sigmay[2]
    sigma_AKT_PP__AKT_t__0_1__model1_data4 = noiseParameter1_AKT_PP__AKT_t__0_1__model1_data4;  // sigmay[3]
    sigma_AKT_PP__AKT_t__0_36__model1_data1 = noiseParameter1_AKT_PP__AKT_t__0_36__model1_data1;  // sigmay[4]
    sigma_AKT_PP__AKT_t__0_69__model1_data1 = noiseParameter1_AKT_PP__AKT_t__0_69__model1_data1;  // sigmay[5]
    sigma_AKT_PP__AKT_t__0_98__model1_data3 = noiseParameter1_AKT_PP__AKT_t__0_98__model1_data3;  // sigmay[6]
    sigma_AKT_PP__AKT_t__1_5__model1_data1 = noiseParameter1_AKT_PP__AKT_t__1_5__model1_data1;  // sigmay[7]
    sigma_AKT_PP__AKT_t__1_6__model1_data4 = noiseParameter1_AKT_PP__AKT_t__1_6__model1_data4;  // sigmay[8]
    sigma_AKT_PP__AKT_t__1_8__model1_data2 = noiseParameter1_AKT_PP__AKT_t__1_8__model1_data2;  // sigmay[9]
    sigma_AKT_PP__AKT_t__2_1__model1_data1 = noiseParameter1_AKT_PP__AKT_t__2_1__model1_data1;  // sigmay[10]
    sigma_AKT_PP__AKT_t__2_3__model1_data1 = noiseParameter1_AKT_PP__AKT_t__2_3__model1_data1;  // sigmay[11]
    sigma_AKT_PP__AKT_t__2_3__model1_data2 = noiseParameter1_AKT_PP__AKT_t__2_3__model1_data2;  // sigmay[12]
    sigma_AKT_PP__AKT_t__2_6__model1_data2 = noiseParameter1_AKT_PP__AKT_t__2_6__model1_data2;  // sigmay[13]
    sigma_AKT_PP__AKT_t__2_8__model1_data2 = noiseParameter1_AKT_PP__AKT_t__2_8__model1_data2;  // sigmay[14]
    sigma_AKT_PP__AKT_t__2_9__model1_data2 = noiseParameter1_AKT_PP__AKT_t__2_9__model1_data2;  // sigmay[15]
    sigma_AKT_PP__AKT_t__3_1__model1_data3 = noiseParameter1_AKT_PP__AKT_t__3_1__model1_data3;  // sigmay[16]
    sigma_AKT_PP__AKT_t__3_2__model1_data2 = noiseParameter1_AKT_PP__AKT_t__3_2__model1_data2;  // sigmay[17]
    sigma_AKT_PP__AKT_t__3_5__model1_data1 = noiseParameter1_AKT_PP__AKT_t__3_5__model1_data1;  // sigmay[18]
    sigma_AKT_PP__AKT_t__3_6__model1_data4 = noiseParameter1_AKT_PP__AKT_t__3_6__model1_data4;  // sigmay[19]
    sigma_AKT_PP__AKT_t__4_5__model1_data3 = noiseParameter1_AKT_PP__AKT_t__4_5__model1_data3;  // sigmay[20]
    sigma_AKT_PP__AKT_t__4_6__model1_data2 = noiseParameter1_AKT_PP__AKT_t__4_6__model1_data2;  // sigmay[21]
    sigma_AKT_PP__AKT_t__4_6__model1_data4 = noiseParameter1_AKT_PP__AKT_t__4_6__model1_data4;  // sigmay[22]
    sigma_AKT_PP__AKT_t__4_7__model1_data2 = noiseParameter1_AKT_PP__AKT_t__4_7__model1_data2;  // sigmay[23]
    sigma_AKT_PP__AKT_t__5_1__model1_data3 = noiseParameter1_AKT_PP__AKT_t__5_1__model1_data3;  // sigmay[24]
    sigma_AKT_PP__AKT_t__5_6__model1_data3 = noiseParameter1_AKT_PP__AKT_t__5_6__model1_data3;  // sigmay[25]
    sigma_AKT_PP__AKT_t__7_4__model1_data1 = noiseParameter1_AKT_PP__AKT_t__7_4__model1_data1;  // sigmay[26]
    sigma_AKT_PP__AKT_t__7_7__model1_data4 = noiseParameter1_AKT_PP__AKT_t__7_7__model1_data4;  // sigmay[27]
    sigma_AKT_PP__AKT_t__9_7__model1_data3 = noiseParameter1_AKT_PP__AKT_t__9_7__model1_data3;  // sigmay[28]
    sigma_AKT_PP__AKT_t__11_1__model1_data3 = noiseParameter1_AKT_PP__AKT_t__11_1__model1_data3;  // sigmay[29]
    sigma_AKT_PP__AKT_t__12_1__model1_data4 = noiseParameter1_AKT_PP__AKT_t__12_1__model1_data4;  // sigmay[30]
    sigma_AKT_PP__AKT_t__13_1__model1_data4 = noiseParameter1_AKT_PP__AKT_t__13_1__model1_data4;  // sigmay[31]
    sigma_AKT_PP__AKT_t__14_1__model1_data3 = noiseParameter1_AKT_PP__AKT_t__14_1__model1_data3;  // sigmay[32]
    sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data1;  // sigmay[33]
    sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data2 = noiseParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data2;  // sigmay[34]
    sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data3 = noiseParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data3;  // sigmay[35]
    sigma_ERB_B1_P_tot__EGFR_t__0_1__model1_data4 = noiseParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data4;  // sigmay[36]
    sigma_ERB_B1_P_tot__EGFR_t__2_4__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__2_4__model1_data1;  // sigmay[37]
    sigma_ERB_B1_P_tot__EGFR_t__2_5__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__2_5__model1_data1;  // sigmay[38]
    sigma_ERB_B1_P_tot__EGFR_t__4_7__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__4_7__model1_data1;  // sigmay[39]
    sigma_ERB_B1_P_tot__EGFR_t__6_8__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__6_8__model1_data1;  // sigmay[40]
    sigma_ERB_B1_P_tot__EGFR_t__7_4__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__7_4__model1_data1;  // sigmay[41]
    sigma_ERB_B1_P_tot__EGFR_t__8_9__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__8_9__model1_data1;  // sigmay[42]
    sigma_ERB_B1_P_tot__EGFR_t__13_1__model1_data1 = noiseParameter1_ERB_B1_P_tot__EGFR_t__13_1__model1_data1;  // sigmay[43]
    sigma_ERK_PP__ERK_t__0_1__model1_data1 = noiseParameter1_ERK_PP__ERK_t__0_1__model1_data1;  // sigmay[44]
    sigma_ERK_PP__ERK_t__0_1__model1_data2 = noiseParameter1_ERK_PP__ERK_t__0_1__model1_data2;  // sigmay[45]
    sigma_ERK_PP__ERK_t__0_1__model1_data3 = noiseParameter1_ERK_PP__ERK_t__0_1__model1_data3;  // sigmay[46]
    sigma_ERK_PP__ERK_t__0_1__model1_data4 = noiseParameter1_ERK_PP__ERK_t__0_1__model1_data4;  // sigmay[47]
    sigma_ERK_PP__ERK_t__0_46__model1_data1 = noiseParameter1_ERK_PP__ERK_t__0_46__model1_data1;  // sigmay[48]
    sigma_ERK_PP__ERK_t__0_94__model1_data2 = noiseParameter1_ERK_PP__ERK_t__0_94__model1_data2;  // sigmay[49]
    sigma_ERK_PP__ERK_t__0_94__model1_data3 = noiseParameter1_ERK_PP__ERK_t__0_94__model1_data3;  // sigmay[50]
    sigma_ERK_PP__ERK_t__1_04__model1_data4 = noiseParameter1_ERK_PP__ERK_t__1_04__model1_data4;  // sigmay[51]
    sigma_ERK_PP__ERK_t__1_1__model1_data4 = noiseParameter1_ERK_PP__ERK_t__1_1__model1_data4;  // sigmay[52]
    sigma_ERK_PP__ERK_t__1_2__model1_data1 = noiseParameter1_ERK_PP__ERK_t__1_2__model1_data1;  // sigmay[53]
    sigma_ERK_PP__ERK_t__1_2__model1_data2 = noiseParameter1_ERK_PP__ERK_t__1_2__model1_data2;  // sigmay[54]
    sigma_ERK_PP__ERK_t__1_2__model1_data3 = noiseParameter1_ERK_PP__ERK_t__1_2__model1_data3;  // sigmay[55]
    sigma_ERK_PP__ERK_t__1_5__model1_data1 = noiseParameter1_ERK_PP__ERK_t__1_5__model1_data1;  // sigmay[56]
    sigma_ERK_PP__ERK_t__1_5__model1_data4 = noiseParameter1_ERK_PP__ERK_t__1_5__model1_data4;  // sigmay[57]
    sigma_ERK_PP__ERK_t__1_7__model1_data4 = noiseParameter1_ERK_PP__ERK_t__1_7__model1_data4;  // sigmay[58]
    sigma_ERK_PP__ERK_t__1_9__model1_data1 = noiseParameter1_ERK_PP__ERK_t__1_9__model1_data1;  // sigmay[59]
    sigma_ERK_PP__ERK_t__1_9__model1_data2 = noiseParameter1_ERK_PP__ERK_t__1_9__model1_data2;  // sigmay[60]
    sigma_ERK_PP__ERK_t__1_9__model1_data3 = noiseParameter1_ERK_PP__ERK_t__1_9__model1_data3;  // sigmay[61]
    sigma_ERK_PP__ERK_t__2_2__model1_data4 = noiseParameter1_ERK_PP__ERK_t__2_2__model1_data4;  // sigmay[62]
    sigma_ERK_PP__ERK_t__2_3__model1_data2 = noiseParameter1_ERK_PP__ERK_t__2_3__model1_data2;  // sigmay[63]
    sigma_ERK_PP__ERK_t__2_3__model1_data3 = noiseParameter1_ERK_PP__ERK_t__2_3__model1_data3;  // sigmay[64]
    sigma_ERK_PP__ERK_t__2_7__model1_data2 = noiseParameter1_ERK_PP__ERK_t__2_7__model1_data2;  // sigmay[65]
    sigma_ERK_PP__ERK_t__2_7__model1_data3 = noiseParameter1_ERK_PP__ERK_t__2_7__model1_data3;  // sigmay[66]
    sigma_ERK_PP__ERK_t__2_8__model1_data1 = noiseParameter1_ERK_PP__ERK_t__2_8__model1_data1;  // sigmay[67]
    sigma_ERK_PP__ERK_t__2_9__model1_data1 = noiseParameter1_ERK_PP__ERK_t__2_9__model1_data1;  // sigmay[68]
    sigma_ERK_PP__ERK_t__3_0__model1_data4 = noiseParameter1_ERK_PP__ERK_t__3_0__model1_data4;  // sigmay[69]
    sigma_ERK_PP__ERK_t__3_1__model1_data1 = noiseParameter1_ERK_PP__ERK_t__3_1__model1_data1;  // sigmay[70]
    sigma_ERK_PP__ERK_t__3_4__model1_data2 = noiseParameter1_ERK_PP__ERK_t__3_4__model1_data2;  // sigmay[71]
    sigma_ERK_PP__ERK_t__3_4__model1_data3 = noiseParameter1_ERK_PP__ERK_t__3_4__model1_data3;  // sigmay[72]
    sigma_ERK_PP__ERK_t__3_7__model1_data1 = noiseParameter1_ERK_PP__ERK_t__3_7__model1_data1;  // sigmay[73]
    sigma_ERK_PP__ERK_t__4_0__model1_data2 = noiseParameter1_ERK_PP__ERK_t__4_0__model1_data2;  // sigmay[74]
    sigma_ERK_PP__ERK_t__4_0__model1_data3 = noiseParameter1_ERK_PP__ERK_t__4_0__model1_data3;  // sigmay[75]
    sigma_ERK_PP__ERK_t__4_3__model1_data2 = noiseParameter1_ERK_PP__ERK_t__4_3__model1_data2;  // sigmay[76]
    sigma_ERK_PP__ERK_t__4_3__model1_data3 = noiseParameter1_ERK_PP__ERK_t__4_3__model1_data3;  // sigmay[77]
    sigma_ERK_PP__ERK_t__4_3__model1_data4 = noiseParameter1_ERK_PP__ERK_t__4_3__model1_data4;  // sigmay[78]
    sigma_ERK_PP__ERK_t__4_7__model1_data4 = noiseParameter1_ERK_PP__ERK_t__4_7__model1_data4;  // sigmay[79]
    sigma_ERK_PP__ERK_t__6_8__model1_data2 = noiseParameter1_ERK_PP__ERK_t__6_8__model1_data2;  // sigmay[80]
    sigma_ERK_PP__ERK_t__6_8__model1_data3 = noiseParameter1_ERK_PP__ERK_t__6_8__model1_data3;  // sigmay[81]
    sigma_ERK_PP__ERK_t__8_2__model1_data4 = noiseParameter1_ERK_PP__ERK_t__8_2__model1_data4;  // sigmay[82]
}

} // namespace model_Chen_MSB2009
} // namespace amici
