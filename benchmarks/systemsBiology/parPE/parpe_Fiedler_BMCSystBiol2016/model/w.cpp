#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void w_Fiedler_BMCSystBiol2016(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    k1max = k10 - k11*(-1 + std::exp(-t/tau1))*std::exp(-t/tau2);  // w[0]
    flux_v1_v_0 = 1.0*K_1*RAF*k1max/(K_1 + pERK);  // w[1]
    flux_v2_v_1 = 1.0*k2*pRAF;  // w[2]
    flux_v3_v_2 = 1.0*K_2*MEK*k3*pRAF/(K_2 + Sorafenib);  // w[3]
    flux_v4_v_3 = 1.0*k4*pMEK;  // w[4]
    flux_v5_v_4 = 1.0*ERK*K_3*k5*pMEK/(K_3 + UO126);  // w[5]
    flux_v6_v_5 = 1.0*k6*pERK;  // w[6]
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
