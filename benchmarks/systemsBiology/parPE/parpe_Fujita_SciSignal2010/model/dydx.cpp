#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Fujita_SciSignal2010 {

void dydx_Fujita_SciSignal2010(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[87] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1;
    dydx[88] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3;
    dydx[89] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0;
    dydx[90] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0;
    dydx[91] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0;
    dydx[92] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0;
    dydx[93] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0;
    dydx[94] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0;
    dydx[95] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0;
    dydx[96] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0;
    dydx[97] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0;
    dydx[98] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0;
    dydx[99] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0;
    dydx[100] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0;
    dydx[101] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0;
    dydx[102] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0;
    dydx[103] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0;
    dydx[104] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0;
    dydx[105] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0;
    dydx[106] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0;
    dydx[107] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0;
    dydx[108] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0;
    dydx[109] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0;
    dydx[154] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1;
    dydx[155] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3;
    dydx[156] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0;
    dydx[157] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0;
    dydx[158] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0;
    dydx[159] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0;
    dydx[160] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0;
    dydx[161] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0;
    dydx[162] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0;
    dydx[163] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0;
    dydx[164] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0;
    dydx[165] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0;
    dydx[166] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0;
    dydx[167] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0;
    dydx[168] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0;
    dydx[169] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0;
    dydx[170] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0;
    dydx[171] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0;
    dydx[172] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0;
    dydx[173] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0;
    dydx[174] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0;
    dydx[175] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0;
    dydx[176] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0;
    dydx[268] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1;
    dydx[269] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3;
    dydx[270] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0;
    dydx[271] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0;
    dydx[272] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0;
    dydx[273] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0;
    dydx[274] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3;
    dydx[275] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0;
    dydx[276] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0;
    dydx[277] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0;
    dydx[278] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0;
    dydx[279] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0;
    dydx[280] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3;
    dydx[281] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0;
    dydx[282] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0;
    dydx[283] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0;
    dydx[284] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0;
    dydx[285] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0;
    dydx[286] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0;
    dydx[287] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0;
    dydx[402] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1;
    dydx[403] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3;
    dydx[404] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0;
    dydx[405] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0;
    dydx[406] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0;
    dydx[407] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0;
    dydx[408] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3;
    dydx[409] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0;
    dydx[410] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0;
    dydx[411] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0;
    dydx[412] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0;
    dydx[413] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0;
    dydx[414] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3;
    dydx[415] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0;
    dydx[416] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0;
    dydx[417] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0;
    dydx[418] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0;
    dydx[419] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0;
    dydx[420] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0;
    dydx[421] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0;
    dydx[512] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1;
    dydx[513] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3;
    dydx[514] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0;
    dydx[515] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0;
    dydx[516] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0;
    dydx[517] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0;
    dydx[518] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1;
    dydx[519] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3;
    dydx[520] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0;
    dydx[521] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0;
    dydx[522] = observableParameter1_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3;
    dydx[523] = observableParameter1_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0;
    dydx[524] = observableParameter1_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3;
    dydx[525] = observableParameter1_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0;
    dydx[526] = observableParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1;
    dydx[527] = observableParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0;
    dydx[528] = observableParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0;
    dydx[529] = observableParameter1_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0;
    dydx[530] = observableParameter1_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0;
    dydx[531] = observableParameter1_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0;
    dydx[532] = observableParameter1_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0;
    dydx[533] = observableParameter1_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0;
    dydx[534] = observableParameter1_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0;
    dydx[535] = observableParameter1_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3;
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
