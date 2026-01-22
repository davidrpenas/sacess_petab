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

namespace amici {
namespace model_Fujita_SciSignal2010 {

void y_Fujita_SciSignal2010(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_1*(pAkt + pAkt_S6);
    y[1] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_00_3*(pAkt + pAkt_S6);
    y[2] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_01_0*(pAkt + pAkt_S6);
    y[3] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_03_0*(pAkt + pAkt_S6);
    y[4] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_10_0*(pAkt + pAkt_S6);
    y[5] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_01__condition_step_30_0*(pAkt + pAkt_S6);
    y[6] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_00_3*(pAkt + pAkt_S6);
    y[7] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_01_0*(pAkt + pAkt_S6);
    y[8] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_03_0*(pAkt + pAkt_S6);
    y[9] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_10_0*(pAkt + pAkt_S6);
    y[10] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_02__condition_step_30_0*(pAkt + pAkt_S6);
    y[11] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_025__condition_step_10_0*(pAkt + pAkt_S6);
    y[12] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_00_3*(pAkt + pAkt_S6);
    y[13] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_01_0*(pAkt + pAkt_S6);
    y[14] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_03__condition_step_03_0*(pAkt + pAkt_S6);
    y[15] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_03_0*(pAkt + pAkt_S6);
    y[16] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_05__condition_step_30_0*(pAkt + pAkt_S6);
    y[17] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_07__condition_step_01_0*(pAkt + pAkt_S6);
    y[18] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_01_0*(pAkt + pAkt_S6);
    y[19] = observableParameter1_pAkt_tot__scaling_pAkt_tot__0_1__condition_step_10_0*(pAkt + pAkt_S6);
    y[20] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_1*(pEGFR + pEGFR_Akt);
    y[21] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_00_3*(pEGFR + pEGFR_Akt);
    y[22] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_01_0*(pEGFR + pEGFR_Akt);
    y[23] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_03_0*(pEGFR + pEGFR_Akt);
    y[24] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_10_0*(pEGFR + pEGFR_Akt);
    y[25] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_01__condition_step_30_0*(pEGFR + pEGFR_Akt);
    y[26] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_01_0*(pEGFR + pEGFR_Akt);
    y[27] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_03_0*(pEGFR + pEGFR_Akt);
    y[28] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_10_0*(pEGFR + pEGFR_Akt);
    y[29] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_02__condition_step_30_0*(pEGFR + pEGFR_Akt);
    y[30] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_025__condition_step_01_0*(pEGFR + pEGFR_Akt);
    y[31] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_03__condition_step_03_0*(pEGFR + pEGFR_Akt);
    y[32] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_01_0*(pEGFR + pEGFR_Akt);
    y[33] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_10_0*(pEGFR + pEGFR_Akt);
    y[34] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_05__condition_step_30_0*(pEGFR + pEGFR_Akt);
    y[35] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_01_0*(pEGFR + pEGFR_Akt);
    y[36] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_065__condition_step_03_0*(pEGFR + pEGFR_Akt);
    y[37] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_1__condition_step_30_0*(pEGFR + pEGFR_Akt);
    y[38] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_12__condition_step_03_0*(pEGFR + pEGFR_Akt);
    y[39] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_03_0*(pEGFR + pEGFR_Akt);
    y[40] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_15__condition_step_10_0*(pEGFR + pEGFR_Akt);
    y[41] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_2__condition_step_30_0*(pEGFR + pEGFR_Akt);
    y[42] = observableParameter1_pEGFR_tot__scaling_pEGFR_tot__0_3__condition_step_10_0*(pEGFR + pEGFR_Akt);
    y[43] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_1*pS6;
    y[44] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_00_3*pS6;
    y[45] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_01_0*pS6;
    y[46] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_03_0*pS6;
    y[47] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_10_0*pS6;
    y[48] = observableParameter1_pS6_tot__scaling_pS6_tot__0_01__condition_step_30_0*pS6;
    y[49] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_1*pS6;
    y[50] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_00_3*pS6;
    y[51] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_10_0*pS6;
    y[52] = observableParameter1_pS6_tot__scaling_pS6_tot__0_02__condition_step_30_0*pS6;
    y[53] = observableParameter1_pS6_tot__scaling_pS6_tot__0_025__condition_step_00_3*pS6;
    y[54] = observableParameter1_pS6_tot__scaling_pS6_tot__0_025__condition_step_01_0*pS6;
    y[55] = observableParameter1_pS6_tot__scaling_pS6_tot__0_03__condition_step_00_3*pS6;
    y[56] = observableParameter1_pS6_tot__scaling_pS6_tot__0_03__condition_step_30_0*pS6;
    y[57] = observableParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_00_1*pS6;
    y[58] = observableParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_03_0*pS6;
    y[59] = observableParameter1_pS6_tot__scaling_pS6_tot__0_035__condition_step_30_0*pS6;
    y[60] = observableParameter1_pS6_tot__scaling_pS6_tot__0_04__condition_step_10_0*pS6;
    y[61] = observableParameter1_pS6_tot__scaling_pS6_tot__0_07__condition_step_01_0*pS6;
    y[62] = observableParameter1_pS6_tot__scaling_pS6_tot__0_07__condition_step_30_0*pS6;
    y[63] = observableParameter1_pS6_tot__scaling_pS6_tot__0_075__condition_step_01_0*pS6;
    y[64] = observableParameter1_pS6_tot__scaling_pS6_tot__0_08__condition_step_03_0*pS6;
    y[65] = observableParameter1_pS6_tot__scaling_pS6_tot__0_08__condition_step_10_0*pS6;
    y[66] = observableParameter1_pS6_tot__scaling_pS6_tot__0_1__condition_step_00_3*pS6;
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
