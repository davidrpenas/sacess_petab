#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 7> dxdotdw_colptrs_Bruno_JExpBot2016_ = {
    0, 3, 5, 8, 11, 13, 16
};

void dxdotdw_colptrs_Bruno_JExpBot2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bruno_JExpBot2016 {

static constexpr std::array<sunindextype, 16> dxdotdw_rowvals_Bruno_JExpBot2016_ = {
    0, 2, 3, 2, 3, 1, 2, 5, 1, 3, 4, 4, 5, 4, 5, 6
};

void dxdotdw_rowvals_Bruno_JExpBot2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Bruno_JExpBot2016_));
}
} // namespace model_Bruno_JExpBot2016
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "dxdotdw.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dxdotdw_Bruno_JExpBot2016(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddbcardt_dflux_v1_ReactionName = -1;  // dxdotdw[0]
    ddb10dt_dflux_v1_ReactionName = 1;  // dxdotdw[1]
    ddbiodt_dflux_v1_ReactionName = 1;  // dxdotdw[2]
    ddb10dt_dflux_v2_ReactionName = -1;  // dxdotdw[3]
    ddbiodt_dflux_v2_ReactionName = 1;  // dxdotdw[4]
    ddbcrydt_dflux_v3_ReactionName = -1;  // dxdotdw[5]
    ddb10dt_dflux_v3_ReactionName = 1;  // dxdotdw[6]
    ddohbiodt_dflux_v3_ReactionName = 1;  // dxdotdw[7]
    ddbcrydt_dflux_v4_ReactionName = -1;  // dxdotdw[8]
    ddbiodt_dflux_v4_ReactionName = 1;  // dxdotdw[9]
    ddohb10dt_dflux_v4_ReactionName = 1;  // dxdotdw[10]
    ddohb10dt_dflux_v5_ReactionName = -1;  // dxdotdw[11]
    ddohbiodt_dflux_v5_ReactionName = 1;  // dxdotdw[12]
    ddohb10dt_dflux_v6_ReactionName = 1;  // dxdotdw[13]
    ddohbiodt_dflux_v6_ReactionName = 1;  // dxdotdw[14]
    ddzeadt_dflux_v6_ReactionName = -1;  // dxdotdw[15]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
