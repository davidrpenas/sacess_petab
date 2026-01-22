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
namespace model_Perelson_Science1996 {

void xdot_Perelson_Science1996(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dTstardt = flux_v1 - flux_v2;  // xdot[0]
    dVdt = -flux_v3 - flux_v4 + flux_v5;  // xdot[1]
    dVindt = -flux_v3;  // xdot[2]
    dVnidt = -flux_v4 + flux_v5;  // xdot[3]
}

} // namespace model_Perelson_Science1996
} // namespace amici
