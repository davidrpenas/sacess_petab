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

void Jy_Fujita_SciSignal2010(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0 + pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0 + pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3 + pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0 + pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0 + pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0 + pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0, 2)) + 0.5*std::pow(-mpAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0 + pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0, 2)/std::pow(sigma_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0, 2);
            break;
        case 23:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0, 2);
            break;
        case 24:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0, 2);
            break;
        case 25:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0, 2);
            break;
        case 26:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0, 2);
            break;
        case 27:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0, 2);
            break;
        case 28:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0, 2);
            break;
        case 29:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0, 2);
            break;
        case 30:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0, 2);
            break;
        case 31:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0, 2);
            break;
        case 32:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0, 2);
            break;
        case 33:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0, 2);
            break;
        case 34:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0, 2);
            break;
        case 35:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0 + pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0, 2);
            break;
        case 36:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0, 2);
            break;
        case 37:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0, 2);
            break;
        case 38:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0, 2);
            break;
        case 39:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0 + pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0, 2);
            break;
        case 40:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0, 2);
            break;
        case 41:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0 + pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0, 2);
            break;
        case 42:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0, 2)) + 0.5*std::pow(-mpEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0 + pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0, 2)/std::pow(sigma_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0, 2);
            break;
        case 43:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_00_1 + pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1, 2);
            break;
        case 44:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3, 2);
            break;
        case 45:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0, 2);
            break;
        case 46:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_03_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0, 2);
            break;
        case 47:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0, 2);
            break;
        case 48:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_01__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0, 2);
            break;
        case 49:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_00_1 + pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1, 2);
            break;
        case 50:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3, 2);
            break;
        case 51:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0, 2);
            break;
        case 52:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_02__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0, 2);
            break;
        case 53:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_025__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3, 2);
            break;
        case 54:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_025__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0, 2);
            break;
        case 55:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_03__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3, 2);
            break;
        case 56:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_03__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0, 2);
            break;
        case 57:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_035__condition_step_00_1 + pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1, 2);
            break;
        case 58:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_035__condition_step_03_0 + pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0, 2);
            break;
        case 59:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_035__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0, 2);
            break;
        case 60:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_04__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0, 2);
            break;
        case 61:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_07__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0, 2);
            break;
        case 62:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_07__condition_step_30_0 + pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0, 2);
            break;
        case 63:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_075__condition_step_01_0 + pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0, 2);
            break;
        case 64:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_08__condition_step_03_0 + pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0, 2);
            break;
        case 65:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_08__condition_step_10_0 + pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0, 2);
            break;
        case 66:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3, 2)) + 0.5*std::pow(-mpS6_tot__scaling_pS6_tot__0_1__condition_step_00_3 + pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3, 2)/std::pow(sigma_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3, 2);
            break;
    }
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
