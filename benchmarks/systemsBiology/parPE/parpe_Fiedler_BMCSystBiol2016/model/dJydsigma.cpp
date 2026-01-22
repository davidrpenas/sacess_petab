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
namespace model_Fiedler_BMCSystBiol2016 {

void dJydsigma_Fiedler_BMCSystBiol2016(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1 - 1.0*std::pow(-mpErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1 + pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1, 2)/std::pow(sigma_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2 - 1.0*std::pow(-mpErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2 + pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2, 2)/std::pow(sigma_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3 - 1.0*std::pow(-mpErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3 + pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3, 2)/std::pow(sigma_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1 - 1.0*std::pow(-mpErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1 + pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1, 2)/std::pow(sigma_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2 - 1.0*std::pow(-mpErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2 + pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2, 2)/std::pow(sigma_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3 - 1.0*std::pow(-mpErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3 + pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3, 2)/std::pow(sigma_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1 - 1.0*std::pow(-mpMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1 + pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1, 2)/std::pow(sigma_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2 - 1.0*std::pow(-mpMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2 + pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2, 2)/std::pow(sigma_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3 - 1.0*std::pow(-mpMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3 + pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3, 2)/std::pow(sigma_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1 - 1.0*std::pow(-mpMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1 + pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1, 2)/std::pow(sigma_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2 - 1.0*std::pow(-mpMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2 + pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2, 2)/std::pow(sigma_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3 - 1.0*std::pow(-mpMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3 + pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3, 2)/std::pow(sigma_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3, 3);
            break;
    }
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
