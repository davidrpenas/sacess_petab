#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<sunindextype, 34> dwdp_colptrs_Crauste_CellSystems2017_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12
};

void dwdp_colptrs_Crauste_CellSystems2017(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Crauste_CellSystems2017_));
}
} // namespace model_Crauste_CellSystems2017
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<sunindextype, 12> dwdp_rowvals_Crauste_CellSystems2017_ = {
    4, 7, 1, 3, 6, 5, 0, 11, 9, 10, 2, 8
};

void dwdp_rowvals_Crauste_CellSystems2017(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Crauste_CellSystems2017_));
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
#include "dwdp.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void dwdp_Crauste_CellSystems2017(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v5_v_4_ddelta_EL = EarlyEffector;  // dwdp[0]
    dflux_v8_v_7_ddelta_LM = LateEffector;  // dwdp[1]
    dflux_v2_v_1_ddelta_NE = Naive*Pathogen;  // dwdp[2]
    dflux_v4_v_3_dmu_EE = std::pow(EarlyEffector, 2);  // dwdp[3]
    dflux_v7_v_6_dmu_LE = EarlyEffector*LateEffector;  // dwdp[4]
    dflux_v6_v_5_dmu_LL = std::pow(LateEffector, 2);  // dwdp[5]
    dflux_v1_v_0_dmu_N = Naive;  // dwdp[6]
    dflux_v12_v_11_dmu_P = Pathogen;  // dwdp[7]
    dflux_v10_v_9_dmu_PE = EarlyEffector*Pathogen;  // dwdp[8]
    dflux_v11_v_10_dmu_PL = LateEffector*Pathogen;  // dwdp[9]
    dflux_v3_v_2_drho_E = EarlyEffector*Pathogen;  // dwdp[10]
    dflux_v9_v_8_drho_P = std::pow(Pathogen, 2);  // dwdp[11]
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
