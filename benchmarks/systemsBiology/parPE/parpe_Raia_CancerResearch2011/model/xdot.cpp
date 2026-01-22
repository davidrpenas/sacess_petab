#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void xdot_Raia_CancerResearch2011(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dRecdt = -0.01*flux_v1_v_0 - 0.01*flux_v2_v_1 + 0.01*flux_v3_v_2;  // xdot[0]
    dRec_idt = 0.01*flux_v2_v_1 - 0.01*flux_v3_v_2;  // xdot[1]
    dIL13_Recdt = 0.01*flux_v1_v_0 - 0.01*flux_v4_v_3;  // xdot[2]
    dp_IL13_Recdt = 0.01*flux_v4_v_3 - 0.01*flux_v7_v_6;  // xdot[3]
    dp_IL13_Rec_idt = 0.01*flux_v7_v_6 - 0.01*flux_v8_v_7;  // xdot[4]
    dJAK2dt = -0.01*flux_v5_v_4 - 0.01*flux_v6_v_5 + 0.01*flux_v9_v_8;  // xdot[5]
    dpJAK2dt = 0.01*flux_v5_v_4 + 0.01*flux_v6_v_5 - 0.01*flux_v9_v_8;  // xdot[6]
    dSTAT5dt = -0.01*flux_v10_v_9 + 0.01*flux_v11_v_10;  // xdot[7]
    dpSTAT5dt = 0.01*flux_v10_v_9 - 0.01*flux_v11_v_10;  // xdot[8]
    dSOCS3mRNAdt = 0.10000000000000001*flux_v13_v_12;  // xdot[9]
    dDecoyRdt = -0.01*flux_v12_v_11;  // xdot[10]
    dIL13_DecoyRdt = 0.01*flux_v12_v_11;  // xdot[11]
    dSOCS3dt = 0.01*flux_v14_v_13 - 0.01*flux_v15_v_14;  // xdot[12]
    dCD274mRNAdt = 0.10000000000000001*flux_v16_v_15;  // xdot[13]
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
