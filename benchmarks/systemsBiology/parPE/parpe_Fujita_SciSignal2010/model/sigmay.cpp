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
namespace model_Fujita_SciSignal2010 {

void sigmay_Fujita_SciSignal2010(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1;  // sigmay[0]
    sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3;  // sigmay[1]
    sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0;  // sigmay[2]
    sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0;  // sigmay[3]
    sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0;  // sigmay[4]
    sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0;  // sigmay[5]
    sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3;  // sigmay[6]
    sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0;  // sigmay[7]
    sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0;  // sigmay[8]
    sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0;  // sigmay[9]
    sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0;  // sigmay[10]
    sigma_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0;  // sigmay[11]
    sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3;  // sigmay[12]
    sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0;  // sigmay[13]
    sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0;  // sigmay[14]
    sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0;  // sigmay[15]
    sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0;  // sigmay[16]
    sigma_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0;  // sigmay[17]
    sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0;  // sigmay[18]
    sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0 = noiseParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0;  // sigmay[19]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1;  // sigmay[20]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3;  // sigmay[21]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0;  // sigmay[22]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0;  // sigmay[23]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0;  // sigmay[24]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0;  // sigmay[25]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0;  // sigmay[26]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0;  // sigmay[27]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0;  // sigmay[28]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0;  // sigmay[29]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0;  // sigmay[30]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0;  // sigmay[31]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0;  // sigmay[32]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0;  // sigmay[33]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0;  // sigmay[34]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0;  // sigmay[35]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0;  // sigmay[36]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0;  // sigmay[37]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0;  // sigmay[38]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0;  // sigmay[39]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0;  // sigmay[40]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0;  // sigmay[41]
    sigma_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0 = noiseParameter1_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0;  // sigmay[42]
    sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1;  // sigmay[43]
    sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3;  // sigmay[44]
    sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0;  // sigmay[45]
    sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0;  // sigmay[46]
    sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0;  // sigmay[47]
    sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0;  // sigmay[48]
    sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1;  // sigmay[49]
    sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3;  // sigmay[50]
    sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0;  // sigmay[51]
    sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0;  // sigmay[52]
    sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3;  // sigmay[53]
    sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0;  // sigmay[54]
    sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3;  // sigmay[55]
    sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0;  // sigmay[56]
    sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1;  // sigmay[57]
    sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0;  // sigmay[58]
    sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0;  // sigmay[59]
    sigma_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0;  // sigmay[60]
    sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0;  // sigmay[61]
    sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0;  // sigmay[62]
    sigma_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0;  // sigmay[63]
    sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0;  // sigmay[64]
    sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0;  // sigmay[65]
    sigma_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3 = noiseParameter1_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3;  // sigmay[66]
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
