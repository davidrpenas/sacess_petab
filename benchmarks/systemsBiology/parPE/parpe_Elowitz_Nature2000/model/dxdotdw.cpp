#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<sunindextype, 17> dxdotdw_colptrs_Elowitz_Nature2000_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16
};

void dxdotdw_colptrs_Elowitz_Nature2000(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Elowitz_Nature2000_));
}
} // namespace model_Elowitz_Nature2000
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<sunindextype, 16> dxdotdw_rowvals_Elowitz_Nature2000_ = {
    4, 5, 6, 7, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 6, 3
};

void dxdotdw_rowvals_Elowitz_Nature2000(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Elowitz_Nature2000_));
}
} // namespace model_Elowitz_Nature2000
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
namespace model_Elowitz_Nature2000 {

void dxdotdw_Elowitz_Nature2000(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddX_mRNAdt_dflux_v1_v_0 = -1;  // dxdotdw[0]
    ddY_mRNAdt_dflux_v2_v_1 = -1;  // dxdotdw[1]
    ddZ_mRNAdt_dflux_v3_v_2 = -1;  // dxdotdw[2]
    ddGFP_mRNAdt_dflux_v4_v_3 = -1;  // dxdotdw[3]
    ddX_proteindt_dflux_v5_v_4 = 1;  // dxdotdw[4]
    ddY_proteindt_dflux_v6_v_5 = 1;  // dxdotdw[5]
    ddZ_proteindt_dflux_v7_v_6 = 1;  // dxdotdw[6]
    ddGFPdt_dflux_v8_v_7 = 1;  // dxdotdw[7]
    ddX_proteindt_dflux_v9_v_8 = -1;  // dxdotdw[8]
    ddY_proteindt_dflux_v10_v_9 = -1;  // dxdotdw[9]
    ddZ_proteindt_dflux_v11_v_10 = -1;  // dxdotdw[10]
    ddGFPdt_dflux_v12_v_11 = -1;  // dxdotdw[11]
    ddX_mRNAdt_dflux_v13_v_12 = 1;  // dxdotdw[12]
    ddY_mRNAdt_dflux_v14_v_13 = 1;  // dxdotdw[13]
    ddZ_mRNAdt_dflux_v15_v_14 = 1;  // dxdotdw[14]
    ddGFPdt_dflux_v16_v_15 = 1;  // dxdotdw[15]
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
