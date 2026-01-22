#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void y_Fiedler_BMCSystBiol2016(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = observableParameter1_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1*pERK;
    y[1] = observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2*pERK;
    y[2] = observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3*pERK;
    y[3] = observableParameter1_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1*pERK;
    y[4] = observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2*pERK;
    y[5] = observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3*pERK;
    y[6] = observableParameter1_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1*pMEK;
    y[7] = observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2*pMEK;
    y[8] = observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3*pMEK;
    y[9] = observableParameter1_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1*pMEK;
    y[10] = observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2*pMEK;
    y[11] = observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3*pMEK;
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
