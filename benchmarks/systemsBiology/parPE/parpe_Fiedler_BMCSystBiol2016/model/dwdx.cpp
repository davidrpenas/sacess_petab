#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 7> dwdx_colptrs_Fiedler_BMCSystBiol2016_ = {
    0, 1, 3, 4, 6, 7, 9
};

void dwdx_colptrs_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Fiedler_BMCSystBiol2016_));
}
} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 9> dwdx_rowvals_Fiedler_BMCSystBiol2016_ = {
    1, 2, 3, 3, 4, 5, 5, 1, 6
};

void dwdx_rowvals_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Fiedler_BMCSystBiol2016_));
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
#include "dwdx.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void dwdx_Fiedler_BMCSystBiol2016(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v2_v_1_dpRAF = 1.0*k2;  // dwdx[1]
        dflux_v4_v_3_dpMEK = 1.0*k4;  // dwdx[4]
        dflux_v6_v_5_dpERK = 1.0*k6;  // dwdx[8]
    }

    // dynamic expressions
    dflux_v1_v_0_dRAF = 1.0*K_1*k1max/(K_1 + pERK);  // dwdx[0]
    dflux_v3_v_2_dpRAF = 1.0*K_2*MEK*k3/(K_2 + Sorafenib);  // dwdx[2]
    dflux_v3_v_2_dMEK = 1.0*K_2*k3*pRAF/(K_2 + Sorafenib);  // dwdx[3]
    dflux_v5_v_4_dpMEK = 1.0*ERK*K_3*k5/(K_3 + UO126);  // dwdx[5]
    dflux_v5_v_4_dERK = 1.0*K_3*k5*pMEK/(K_3 + UO126);  // dwdx[6]
    dflux_v1_v_0_dpERK = -1.0*K_1*RAF*k1max/std::pow(K_1 + pERK, 2);  // dwdx[7]
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
