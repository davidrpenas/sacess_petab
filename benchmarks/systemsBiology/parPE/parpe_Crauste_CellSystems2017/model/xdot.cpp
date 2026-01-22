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
namespace model_Crauste_CellSystems2017 {

void xdot_Crauste_CellSystems2017(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dNaivedt = -flux_v1_v_0 - flux_v2_v_1;  // xdot[0]
    dEarlyEffectordt = flux_v2_v_1 + flux_v3_v_2 - flux_v4_v_3 - flux_v5_v_4;  // xdot[1]
    dLateEffectordt = flux_v5_v_4 - flux_v6_v_5 - flux_v7_v_6 - flux_v8_v_7;  // xdot[2]
    dMemorydt = flux_v8_v_7;  // xdot[3]
    dPathogendt = -flux_v10_v_9 - flux_v11_v_10 - flux_v12_v_11 + flux_v9_v_8;  // xdot[4]
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
