#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<sunindextype, 6> dwdx_colptrs_Crauste_CellSystems2017_ = {
    0, 2, 7, 11, 11, 17
};

void dwdx_colptrs_Crauste_CellSystems2017(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Crauste_CellSystems2017_));
}
} // namespace model_Crauste_CellSystems2017
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<sunindextype, 17> dwdx_rowvals_Crauste_CellSystems2017_ = {
    0, 1, 2, 3, 4, 6, 9, 5, 6, 7, 10, 1, 2, 8, 9, 10, 11
};

void dwdx_rowvals_Crauste_CellSystems2017(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Crauste_CellSystems2017_));
}
} // namespace model_Crauste_CellSystems2017
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void dwdx_Crauste_CellSystems2017(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v1_v_0_dNaive = mu_N;  // dwdx[0]
        dflux_v5_v_4_dEarlyEffector = delta_EL;  // dwdx[4]
        dflux_v8_v_7_dLateEffector = delta_LM;  // dwdx[9]
        dflux_v12_v_11_dPathogen = mu_P;  // dwdx[16]
    }

    // dynamic expressions
    dflux_v2_v_1_dNaive = Pathogen*delta_NE;  // dwdx[1]
    dflux_v3_v_2_dEarlyEffector = Pathogen*rho_E;  // dwdx[2]
    dflux_v4_v_3_dEarlyEffector = 2*EarlyEffector*mu_EE;  // dwdx[3]
    dflux_v7_v_6_dEarlyEffector = LateEffector*mu_LE;  // dwdx[5]
    dflux_v10_v_9_dEarlyEffector = Pathogen*mu_PE;  // dwdx[6]
    dflux_v6_v_5_dLateEffector = 2*LateEffector*mu_LL;  // dwdx[7]
    dflux_v7_v_6_dLateEffector = EarlyEffector*mu_LE;  // dwdx[8]
    dflux_v11_v_10_dLateEffector = Pathogen*mu_PL;  // dwdx[10]
    dflux_v2_v_1_dPathogen = Naive*delta_NE;  // dwdx[11]
    dflux_v3_v_2_dPathogen = EarlyEffector*rho_E;  // dwdx[12]
    dflux_v9_v_8_dPathogen = 2*Pathogen*rho_P;  // dwdx[13]
    dflux_v10_v_9_dPathogen = EarlyEffector*mu_PE;  // dwdx[14]
    dflux_v11_v_10_dPathogen = LateEffector*mu_PL;  // dwdx[15]
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
