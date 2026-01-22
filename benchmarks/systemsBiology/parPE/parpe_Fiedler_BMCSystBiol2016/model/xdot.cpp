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
namespace model_Fiedler_BMCSystBiol2016 {

void xdot_Fiedler_BMCSystBiol2016(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dRAFdt = -flux_v1_v_0 + flux_v2_v_1;  // xdot[0]
    dpRAFdt = flux_v1_v_0 - flux_v2_v_1;  // xdot[1]
    dMEKdt = -flux_v3_v_2 + flux_v4_v_3;  // xdot[2]
    dpMEKdt = flux_v3_v_2 - flux_v4_v_3;  // xdot[3]
    dERKdt = -flux_v5_v_4 + flux_v6_v_5;  // xdot[4]
    dpERKdt = flux_v5_v_4 - flux_v6_v_5;  // xdot[5]
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
