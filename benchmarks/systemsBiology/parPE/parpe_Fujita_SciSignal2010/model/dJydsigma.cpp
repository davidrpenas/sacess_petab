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
namespace model_Fujita_SciSignal2010 {

void dJydsigma_Fujita_SciSignal2010(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3 + pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0 + pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0 - 1.0*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0 - 1.0*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_00_1 + pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_03_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_00_1 + pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_025__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_025__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0, 3);
            break;
        case 55:
            dJydsigma[55] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_03__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3, 3);
            break;
        case 56:
            dJydsigma[56] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_03__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0, 3);
            break;
        case 57:
            dJydsigma[57] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_035__condition_step_00_1 + pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1, 3);
            break;
        case 58:
            dJydsigma[58] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_035__condition_step_03_0 + pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0, 3);
            break;
        case 59:
            dJydsigma[59] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_035__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0, 3);
            break;
        case 60:
            dJydsigma[60] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_04__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0, 3);
            break;
        case 61:
            dJydsigma[61] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_07__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0, 3);
            break;
        case 62:
            dJydsigma[62] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_07__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0, 3);
            break;
        case 63:
            dJydsigma[63] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_075__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0, 3);
            break;
        case 64:
            dJydsigma[64] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_08__condition_step_03_0 + pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0, 3);
            break;
        case 65:
            dJydsigma[65] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_08__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0, 3);
            break;
        case 66:
            dJydsigma[66] = 1.0/sigma_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3 - 1.0*std::pow(-mpS6_tot__scaling_pS6_tot__0_1__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3, 3);
            break;
    }
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
