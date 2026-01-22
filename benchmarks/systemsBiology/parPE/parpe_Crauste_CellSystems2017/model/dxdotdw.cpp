#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<sunindextype, 13> dxdotdw_colptrs_Crauste_CellSystems2017_ = {
    0, 1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15
};

void dxdotdw_colptrs_Crauste_CellSystems2017(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Crauste_CellSystems2017_));
}
} // namespace model_Crauste_CellSystems2017
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<sunindextype, 15> dxdotdw_rowvals_Crauste_CellSystems2017_ = {
    0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4
};

void dxdotdw_rowvals_Crauste_CellSystems2017(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Crauste_CellSystems2017_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void dxdotdw_Crauste_CellSystems2017(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddNaivedt_dflux_v1_v_0 = -1;  // dxdotdw[0]
    ddNaivedt_dflux_v2_v_1 = -1;  // dxdotdw[1]
    ddEarlyEffectordt_dflux_v2_v_1 = 1;  // dxdotdw[2]
    ddEarlyEffectordt_dflux_v3_v_2 = 1;  // dxdotdw[3]
    ddEarlyEffectordt_dflux_v4_v_3 = -1;  // dxdotdw[4]
    ddEarlyEffectordt_dflux_v5_v_4 = -1;  // dxdotdw[5]
    ddLateEffectordt_dflux_v5_v_4 = 1;  // dxdotdw[6]
    ddLateEffectordt_dflux_v6_v_5 = -1;  // dxdotdw[7]
    ddLateEffectordt_dflux_v7_v_6 = -1;  // dxdotdw[8]
    ddLateEffectordt_dflux_v8_v_7 = -1;  // dxdotdw[9]
    ddMemorydt_dflux_v8_v_7 = 1;  // dxdotdw[10]
    ddPathogendt_dflux_v9_v_8 = 1;  // dxdotdw[11]
    ddPathogendt_dflux_v10_v_9 = -1;  // dxdotdw[12]
    ddPathogendt_dflux_v11_v_10 = -1;  // dxdotdw[13]
    ddPathogendt_dflux_v12_v_11 = -1;  // dxdotdw[14]
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
