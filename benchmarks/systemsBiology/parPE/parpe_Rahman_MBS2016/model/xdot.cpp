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
namespace model_Rahman_MBS2016 {

void xdot_Rahman_MBS2016(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dsusceptibledt = flux_v1 - flux_v10 - flux_v2;  // xdot[0]
    dinfected_normaldt = -flux_v11 + flux_v2 - flux_v3 - flux_v7;  // xdot[1]
    dinfected_moderatedt = -flux_v12 + flux_v3 - flux_v4 - flux_v8;  // xdot[2]
    dinfected_weakdt = -flux_v13 + flux_v4 - flux_v9;  // xdot[3]
    dtreated_normaldt = -flux_v14 + flux_v6 + flux_v7;  // xdot[4]
    dtreated_moderatedt = -flux_v15 + flux_v5 - flux_v6 + flux_v8;  // xdot[5]
    dtreated_weakdt = -flux_v16 - flux_v5 + flux_v9;  // xdot[6]
}

} // namespace model_Rahman_MBS2016
} // namespace amici
