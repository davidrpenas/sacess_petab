#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void dsigmaydp_Fiedler_BMCSystBiol2016(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 24:
            dsigmaydp[0] = noiseParameter2_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1;
            break;
        case 25:
            dsigmaydp[0] = noiseParameter1_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1;
            break;
        case 26:
            dsigmaydp[1] = noiseParameter2_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2;
            break;
        case 27:
            dsigmaydp[1] = noiseParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2;
            break;
        case 28:
            dsigmaydp[2] = noiseParameter2_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3;
            break;
        case 29:
            dsigmaydp[2] = noiseParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3;
            break;
        case 30:
            dsigmaydp[3] = noiseParameter2_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1;
            break;
        case 31:
            dsigmaydp[3] = noiseParameter1_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1;
            break;
        case 32:
            dsigmaydp[4] = noiseParameter2_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2;
            break;
        case 33:
            dsigmaydp[4] = noiseParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2;
            break;
        case 34:
            dsigmaydp[5] = noiseParameter2_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3;
            break;
        case 35:
            dsigmaydp[5] = noiseParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3;
            break;
        case 36:
            dsigmaydp[6] = noiseParameter2_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1;
            break;
        case 37:
            dsigmaydp[6] = noiseParameter1_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1;
            break;
        case 38:
            dsigmaydp[7] = noiseParameter2_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2;
            break;
        case 39:
            dsigmaydp[7] = noiseParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2;
            break;
        case 40:
            dsigmaydp[8] = noiseParameter2_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3;
            break;
        case 41:
            dsigmaydp[8] = noiseParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3;
            break;
        case 42:
            dsigmaydp[9] = noiseParameter2_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1;
            break;
        case 43:
            dsigmaydp[9] = noiseParameter1_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1;
            break;
        case 44:
            dsigmaydp[10] = noiseParameter2_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2;
            break;
        case 45:
            dsigmaydp[10] = noiseParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2;
            break;
        case 46:
            dsigmaydp[11] = noiseParameter2_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3;
            break;
        case 47:
            dsigmaydp[11] = noiseParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3;
            break;
    }
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
