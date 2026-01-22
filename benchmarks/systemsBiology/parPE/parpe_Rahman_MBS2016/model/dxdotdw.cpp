#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 22> dxdotdw_colptrs_Rahman_MBS2016_ = {
    0, 0, 0, 0, 0, 0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 18, 19, 20, 21, 22, 23, 24
};

void dxdotdw_colptrs_Rahman_MBS2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Rahman_MBS2016_));
}
} // namespace model_Rahman_MBS2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 24> dxdotdw_rowvals_Rahman_MBS2016_ = {
    0, 0, 1, 1, 2, 2, 3, 5, 6, 4, 5, 1, 4, 2, 5, 3, 6, 0, 1, 2, 3, 4, 5, 6
};

void dxdotdw_rowvals_Rahman_MBS2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Rahman_MBS2016_));
}
} // namespace model_Rahman_MBS2016
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
namespace model_Rahman_MBS2016 {

void dxdotdw_Rahman_MBS2016(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddsusceptibledt_dflux_v1 = 1;  // dxdotdw[0]
    ddsusceptibledt_dflux_v2 = -1;  // dxdotdw[1]
    ddinfected_normaldt_dflux_v2 = 1;  // dxdotdw[2]
    ddinfected_normaldt_dflux_v3 = -1;  // dxdotdw[3]
    ddinfected_moderatedt_dflux_v3 = 1;  // dxdotdw[4]
    ddinfected_moderatedt_dflux_v4 = -1;  // dxdotdw[5]
    ddinfected_weakdt_dflux_v4 = 1;  // dxdotdw[6]
    ddtreated_moderatedt_dflux_v5 = 1;  // dxdotdw[7]
    ddtreated_weakdt_dflux_v5 = -1;  // dxdotdw[8]
    ddtreated_normaldt_dflux_v6 = 1;  // dxdotdw[9]
    ddtreated_moderatedt_dflux_v6 = -1;  // dxdotdw[10]
    ddinfected_normaldt_dflux_v7 = -1;  // dxdotdw[11]
    ddtreated_normaldt_dflux_v7 = 1;  // dxdotdw[12]
    ddinfected_moderatedt_dflux_v8 = -1;  // dxdotdw[13]
    ddtreated_moderatedt_dflux_v8 = 1;  // dxdotdw[14]
    ddinfected_weakdt_dflux_v9 = -1;  // dxdotdw[15]
    ddtreated_weakdt_dflux_v9 = 1;  // dxdotdw[16]
    ddsusceptibledt_dflux_v10 = -1;  // dxdotdw[17]
    ddinfected_normaldt_dflux_v11 = -1;  // dxdotdw[18]
    ddinfected_moderatedt_dflux_v12 = -1;  // dxdotdw[19]
    ddinfected_weakdt_dflux_v13 = -1;  // dxdotdw[20]
    ddtreated_normaldt_dflux_v14 = -1;  // dxdotdw[21]
    ddtreated_moderatedt_dflux_v15 = -1;  // dxdotdw[22]
    ddtreated_weakdt_dflux_v16 = -1;  // dxdotdw[23]
}

} // namespace model_Rahman_MBS2016
} // namespace amici
