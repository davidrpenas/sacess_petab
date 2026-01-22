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
namespace model_Sneyd_PNAS2002 {

void xdot_Sneyd_PNAS2002(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dIPR_Odt = -flux_v_0 + flux_v_1 - flux_v_4 + flux_v_5 - flux_v_8 + flux_v_9;  // xdot[0]
    dIPR_Rdt = flux_v_0 - flux_v_1 - flux_v_2 + flux_v_3;  // xdot[1]
    dIPR_I1dt = flux_v_2 - flux_v_3;  // xdot[2]
    dIPR_Sdt = flux_v_8 - flux_v_9;  // xdot[3]
    dIPR_Adt = flux_v_4 - flux_v_5 - flux_v_6 + flux_v_7;  // xdot[4]
    dIPR_I2dt = flux_v_6 - flux_v_7;  // xdot[5]
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
