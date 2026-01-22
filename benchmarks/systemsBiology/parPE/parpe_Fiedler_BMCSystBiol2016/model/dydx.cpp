#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void dydx_Fiedler_BMCSystBiol2016(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[42] = observableParameter1_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1;
    dydx[43] = observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2;
    dydx[44] = observableParameter1_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3;
    dydx[45] = observableParameter1_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1;
    dydx[46] = observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2;
    dydx[47] = observableParameter1_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3;
    dydx[60] = observableParameter1_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1;
    dydx[61] = observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2;
    dydx[62] = observableParameter1_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3;
    dydx[63] = observableParameter1_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1;
    dydx[64] = observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2;
    dydx[65] = observableParameter1_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3;
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
