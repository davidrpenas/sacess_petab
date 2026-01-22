#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Borghans_BiophysChem1997 {

static constexpr std::array<sunindextype, 9> dxdotdw_colptrs_Borghans_BiophysChem1997_ = {
    0, 1, 3, 5, 7, 8, 9, 10, 11
};

void dxdotdw_colptrs_Borghans_BiophysChem1997(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Borghans_BiophysChem1997_));
}
} // namespace model_Borghans_BiophysChem1997
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Borghans_BiophysChem1997 {

static constexpr std::array<sunindextype, 11> dxdotdw_rowvals_Borghans_BiophysChem1997_ = {
    0, 0, 1, 0, 1, 0, 1, 0, 2, 2, 2
};

void dxdotdw_rowvals_Borghans_BiophysChem1997(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Borghans_BiophysChem1997_));
}
} // namespace model_Borghans_BiophysChem1997
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
namespace model_Borghans_BiophysChem1997 {

void dxdotdw_Borghans_BiophysChem1997(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddZ_statedt_dflux_v1_v_0 = 1;  // dxdotdw[0]
    ddZ_statedt_dflux_v2_v_1 = -1;  // dxdotdw[1]
    ddY_statedt_dflux_v2_v_1 = 1;  // dxdotdw[2]
    ddZ_statedt_dflux_v3_v_2 = 1;  // dxdotdw[3]
    ddY_statedt_dflux_v3_v_2 = -1;  // dxdotdw[4]
    ddZ_statedt_dflux_v4_v_3 = 1;  // dxdotdw[5]
    ddY_statedt_dflux_v4_v_3 = -1;  // dxdotdw[6]
    ddZ_statedt_dflux_v5_v_4 = -1;  // dxdotdw[7]
    ddA_statedt_dflux_v6_v_5 = 1;  // dxdotdw[8]
    ddA_statedt_dflux_v7_v_6 = -1;  // dxdotdw[9]
    ddA_statedt_dflux_v8_v_7 = -1;  // dxdotdw[10]
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
