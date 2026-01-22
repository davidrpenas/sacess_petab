#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 15> dxdotdw_colptrs_Schwen_PONE2014_ = {
    0, 3, 6, 8, 10, 13, 16, 18, 20, 22, 24, 27, 30, 31, 32
};

void dxdotdw_colptrs_Schwen_PONE2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<sunindextype, 32> dxdotdw_rowvals_Schwen_PONE2014_ = {
    0, 1, 3, 0, 2, 4, 0, 10, 0, 10, 0, 1, 3, 0, 2, 4, 3, 5, 4, 6, 3, 5, 4, 6, 1, 5, 9, 2, 6, 9, 7, 8
};

void dxdotdw_rowvals_Schwen_PONE2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Schwen_PONE2014_));
}
} // namespace model_Schwen_PONE2014
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
namespace model_Schwen_PONE2014 {

void dxdotdw_Schwen_PONE2014(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddInsdt_dflux_v1_ka1 = -1;  // dxdotdw[0]
    ddRec1dt_dflux_v1_ka1 = -1;  // dxdotdw[1]
    ddIR1dt_dflux_v1_ka1 = 1;  // dxdotdw[2]
    ddInsdt_dflux_v2_ka2fold = -1;  // dxdotdw[3]
    ddRec2dt_dflux_v2_ka2fold = -1;  // dxdotdw[4]
    ddIR2dt_dflux_v2_ka2fold = 1;  // dxdotdw[5]
    ddInsdt_dflux_v3_v_2 = -1;  // dxdotdw[6]
    ddBoundUnspecdt_dflux_v3_v_2 = 1;  // dxdotdw[7]
    ddInsdt_dflux_v4_v_3 = 1;  // dxdotdw[8]
    ddBoundUnspecdt_dflux_v4_v_3 = -1;  // dxdotdw[9]
    ddInsdt_dflux_v5_kd1 = 1;  // dxdotdw[10]
    ddRec1dt_dflux_v5_kd1 = 1;  // dxdotdw[11]
    ddIR1dt_dflux_v5_kd1 = -1;  // dxdotdw[12]
    ddInsdt_dflux_v6_kd2fold = 1;  // dxdotdw[13]
    ddRec2dt_dflux_v6_kd2fold = 1;  // dxdotdw[14]
    ddIR2dt_dflux_v6_kd2fold = -1;  // dxdotdw[15]
    ddIR1dt_dflux_v7_v_6 = -1;  // dxdotdw[16]
    ddIR1indt_dflux_v7_v_6 = 1;  // dxdotdw[17]
    ddIR2dt_dflux_v8_v_7 = -1;  // dxdotdw[18]
    ddIR2indt_dflux_v8_v_7 = 1;  // dxdotdw[19]
    ddIR1dt_dflux_v9_v_8 = 1;  // dxdotdw[20]
    ddIR1indt_dflux_v9_v_8 = -1;  // dxdotdw[21]
    ddIR2dt_dflux_v10_v_9 = 1;  // dxdotdw[22]
    ddIR2indt_dflux_v10_v_9 = -1;  // dxdotdw[23]
    ddRec1dt_dflux_v11_v_10 = 1;  // dxdotdw[24]
    ddIR1indt_dflux_v11_v_10 = -1;  // dxdotdw[25]
    ddInsulinFragmentsdt_dflux_v11_v_10 = 1;  // dxdotdw[26]
    ddRec2dt_dflux_v12_v_11 = 1;  // dxdotdw[27]
    ddIR2indt_dflux_v12_v_11 = -1;  // dxdotdw[28]
    ddInsulinFragmentsdt_dflux_v12_v_11 = 1;  // dxdotdw[29]
    ddUptake1dt_dflux_v13_v_12 = 1;  // dxdotdw[30]
    ddUptake2dt_dflux_v14_v_13 = 1;  // dxdotdw[31]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
