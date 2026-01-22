#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 8> dxdotdw_colptrs_Fiedler_BMCSystBiol2016_ = {
    0, 0, 2, 4, 6, 8, 10, 12
};

void dxdotdw_colptrs_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Fiedler_BMCSystBiol2016_));
}
} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 12> dxdotdw_rowvals_Fiedler_BMCSystBiol2016_ = {
    0, 1, 0, 1, 2, 3, 2, 3, 4, 5, 4, 5
};

void dxdotdw_rowvals_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Fiedler_BMCSystBiol2016_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void dxdotdw_Fiedler_BMCSystBiol2016(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddRAFdt_dflux_v1_v_0 = -1;  // dxdotdw[0]
    ddpRAFdt_dflux_v1_v_0 = 1;  // dxdotdw[1]
    ddRAFdt_dflux_v2_v_1 = 1;  // dxdotdw[2]
    ddpRAFdt_dflux_v2_v_1 = -1;  // dxdotdw[3]
    ddMEKdt_dflux_v3_v_2 = -1;  // dxdotdw[4]
    ddpMEKdt_dflux_v3_v_2 = 1;  // dxdotdw[5]
    ddMEKdt_dflux_v4_v_3 = 1;  // dxdotdw[6]
    ddpMEKdt_dflux_v4_v_3 = -1;  // dxdotdw[7]
    ddERKdt_dflux_v5_v_4 = -1;  // dxdotdw[8]
    ddpERKdt_dflux_v5_v_4 = 1;  // dxdotdw[9]
    ddERKdt_dflux_v6_v_5 = 1;  // dxdotdw[10]
    ddpERKdt_dflux_v6_v_5 = -1;  // dxdotdw[11]
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
