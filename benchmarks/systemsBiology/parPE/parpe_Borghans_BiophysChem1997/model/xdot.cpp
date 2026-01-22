#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void xdot_Borghans_BiophysChem1997(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dZ_statedt = flux_v1_v_0 - flux_v2_v_1 + flux_v3_v_2 + flux_v4_v_3 - flux_v5_v_4;  // xdot[0]
    dY_statedt = flux_v2_v_1 - flux_v3_v_2 - flux_v4_v_3;  // xdot[1]
    dA_statedt = flux_v6_v_5 - flux_v7_v_6 - flux_v8_v_7;  // xdot[2]
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
