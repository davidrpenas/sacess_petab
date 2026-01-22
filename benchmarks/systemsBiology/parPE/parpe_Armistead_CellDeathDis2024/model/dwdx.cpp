#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<sunindextype, 5> dwdx_colptrs_Armistead_CellDeathDis2024_ = {
    0, 2, 4, 6, 7
};

void dwdx_colptrs_Armistead_CellDeathDis2024(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Armistead_CellDeathDis2024_));
}
} // namespace model_Armistead_CellDeathDis2024
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<sunindextype, 7> dwdx_rowvals_Armistead_CellDeathDis2024_ = {
    2, 4, 3, 5, 6, 7, 1
};

void dwdx_rowvals_Armistead_CellDeathDis2024(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Armistead_CellDeathDis2024_));
}
} // namespace model_Armistead_CellDeathDis2024
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
namespace model_Armistead_CellDeathDis2024 {

void dwdx_Armistead_CellDeathDis2024(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_R3_dCer = k1;  // dwdx[0]
        dflux_R5_dCer = k_d;  // dwdx[1]
        dflux_R4_dSphingo = k2;  // dwdx[2]
        dflux_R6_dSphingo = k3*(-S_on*alpha_hai1a + 1);  // dwdx[3]
        dflux_R7_dS1P = k4;  // dwdx[4]
        dflux_R8_dS1P = k5;  // dwdx[5]
        dflux_R2_dSphinga = k0;  // dwdx[6]
    }
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
