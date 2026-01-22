#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 49> dwdp_colptrs_Fiedler_BMCSystBiol2016_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12
};

void dwdp_colptrs_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Fiedler_BMCSystBiol2016_));
}
} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 12> dwdp_rowvals_Fiedler_BMCSystBiol2016_ = {
    1, 3, 5, 0, 0, 2, 3, 4, 5, 6, 0, 0
};

void dwdp_rowvals_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Fiedler_BMCSystBiol2016_));
}
} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdp.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void dwdp_Fiedler_BMCSystBiol2016(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dk1max_dk10 = 1;  // dwdp[3]
    }

    // dynamic expressions
    dflux_v1_v_0_dK_1 = -1.0*K_1*RAF*k1max/std::pow(K_1 + pERK, 2) + 1.0*RAF*k1max/(K_1 + pERK);  // dwdp[0]
    dflux_v3_v_2_dK_2 = -1.0*K_2*MEK*k3*pRAF/std::pow(K_2 + Sorafenib, 2) + 1.0*MEK*k3*pRAF/(K_2 + Sorafenib);  // dwdp[1]
    dflux_v5_v_4_dK_3 = -1.0*ERK*K_3*k5*pMEK/std::pow(K_3 + UO126, 2) + 1.0*ERK*k5*pMEK/(K_3 + UO126);  // dwdp[2]
    dk1max_dk11 = (1 - std::exp(-t/tau1))*std::exp(-t/tau2);  // dwdp[4]
    dflux_v2_v_1_dk2 = 1.0*pRAF;  // dwdp[5]
    dflux_v3_v_2_dk3 = 1.0*K_2*MEK*pRAF/(K_2 + Sorafenib);  // dwdp[6]
    dflux_v4_v_3_dk4 = 1.0*pMEK;  // dwdp[7]
    dflux_v5_v_4_dk5 = 1.0*ERK*K_3*pMEK/(K_3 + UO126);  // dwdp[8]
    dflux_v6_v_5_dk6 = 1.0*pERK;  // dwdp[9]
    dk1max_dtau1 = -k11*t*std::exp(-t*(1.0/tau2 + 1.0/tau1))/std::pow(tau1, 2);  // dwdp[10]
    dk1max_dtau2 = -k11*t*(-1 + std::exp(-t/tau1))*std::exp(-t/tau2)/std::pow(tau2, 2);  // dwdp[11]
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
