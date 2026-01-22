#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 58> dxdotdw_colptrs_Blasi_CellSystems2016_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64
};

void dxdotdw_colptrs_Blasi_CellSystems2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Blasi_CellSystems2016_));
}
} // namespace model_Blasi_CellSystems2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 64> dxdotdw_rowvals_Blasi_CellSystems2016_ = {
    0, 2, 0, 1, 0, 3, 0, 4, 2, 5, 2, 6, 2, 7, 1, 5, 1, 8, 1, 9, 3, 6, 3, 8, 3, 10, 4, 7, 4, 9, 4, 10, 5, 11, 5, 12, 6, 11, 6, 13, 7, 12, 7, 13, 8, 11, 8, 14, 9, 12, 9, 14, 10, 13, 10, 14, 11, 15, 12, 15, 13, 15, 14, 15
};

void dxdotdw_rowvals_Blasi_CellSystems2016(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Blasi_CellSystems2016_));
}
} // namespace model_Blasi_CellSystems2016
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
namespace model_Blasi_CellSystems2016 {

void dxdotdw_Blasi_CellSystems2016(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddx_0acdt_dflux_re1 = -1;  // dxdotdw[0]
    ddx_k5dt_dflux_re1 = 1;  // dxdotdw[1]
    ddx_0acdt_dflux_re2 = -1;  // dxdotdw[2]
    ddx_k8dt_dflux_re2 = 1;  // dxdotdw[3]
    ddx_0acdt_dflux_re3 = -1;  // dxdotdw[4]
    ddx_k12dt_dflux_re3 = 1;  // dxdotdw[5]
    ddx_0acdt_dflux_re4 = -1;  // dxdotdw[6]
    ddx_k16dt_dflux_re4 = 1;  // dxdotdw[7]
    ddx_k5dt_dflux_re5 = -1;  // dxdotdw[8]
    ddx_k5k8dt_dflux_re5 = 1;  // dxdotdw[9]
    ddx_k5dt_dflux_re6 = -1;  // dxdotdw[10]
    ddx_k5k12dt_dflux_re6 = 1;  // dxdotdw[11]
    ddx_k5dt_dflux_re7 = -1;  // dxdotdw[12]
    ddx_k5k16dt_dflux_re7 = 1;  // dxdotdw[13]
    ddx_k8dt_dflux_re8 = -1;  // dxdotdw[14]
    ddx_k5k8dt_dflux_re8 = 1;  // dxdotdw[15]
    ddx_k8dt_dflux_re9 = -1;  // dxdotdw[16]
    ddx_k8k12dt_dflux_re9 = 1;  // dxdotdw[17]
    ddx_k8dt_dflux_re10 = -1;  // dxdotdw[18]
    ddx_k8k16dt_dflux_re10 = 1;  // dxdotdw[19]
    ddx_k12dt_dflux_re11 = -1;  // dxdotdw[20]
    ddx_k5k12dt_dflux_re11 = 1;  // dxdotdw[21]
    ddx_k12dt_dflux_re12 = -1;  // dxdotdw[22]
    ddx_k8k12dt_dflux_re12 = 1;  // dxdotdw[23]
    ddx_k12dt_dflux_re13 = -1;  // dxdotdw[24]
    ddx_k12k16dt_dflux_re13 = 1;  // dxdotdw[25]
    ddx_k16dt_dflux_re14 = -1;  // dxdotdw[26]
    ddx_k5k16dt_dflux_re14 = 1;  // dxdotdw[27]
    ddx_k16dt_dflux_re15 = -1;  // dxdotdw[28]
    ddx_k8k16dt_dflux_re15 = 1;  // dxdotdw[29]
    ddx_k16dt_dflux_re16 = -1;  // dxdotdw[30]
    ddx_k12k16dt_dflux_re16 = 1;  // dxdotdw[31]
    ddx_k5k8dt_dflux_re17 = -1;  // dxdotdw[32]
    ddx_k5k8k12dt_dflux_re17 = 1;  // dxdotdw[33]
    ddx_k5k8dt_dflux_re18 = -1;  // dxdotdw[34]
    ddx_k5k8k16dt_dflux_re18 = 1;  // dxdotdw[35]
    ddx_k5k12dt_dflux_re19 = -1;  // dxdotdw[36]
    ddx_k5k8k12dt_dflux_re19 = 1;  // dxdotdw[37]
    ddx_k5k12dt_dflux_re20 = -1;  // dxdotdw[38]
    ddx_k5k12k16dt_dflux_re20 = 1;  // dxdotdw[39]
    ddx_k5k16dt_dflux_re21 = -1;  // dxdotdw[40]
    ddx_k5k8k16dt_dflux_re21 = 1;  // dxdotdw[41]
    ddx_k5k16dt_dflux_re22 = -1;  // dxdotdw[42]
    ddx_k5k12k16dt_dflux_re22 = 1;  // dxdotdw[43]
    ddx_k8k12dt_dflux_re23 = -1;  // dxdotdw[44]
    ddx_k5k8k12dt_dflux_re23 = 1;  // dxdotdw[45]
    ddx_k8k12dt_dflux_re24 = -1;  // dxdotdw[46]
    ddx_k8k12k16dt_dflux_re24 = 1;  // dxdotdw[47]
    ddx_k8k16dt_dflux_re25 = -1;  // dxdotdw[48]
    ddx_k5k8k16dt_dflux_re25 = 1;  // dxdotdw[49]
    ddx_k8k16dt_dflux_re26 = -1;  // dxdotdw[50]
    ddx_k8k12k16dt_dflux_re26 = 1;  // dxdotdw[51]
    ddx_k12k16dt_dflux_re27 = -1;  // dxdotdw[52]
    ddx_k5k12k16dt_dflux_re27 = 1;  // dxdotdw[53]
    ddx_k12k16dt_dflux_re28 = -1;  // dxdotdw[54]
    ddx_k8k12k16dt_dflux_re28 = 1;  // dxdotdw[55]
    ddx_k5k8k12dt_dflux_re29 = -1;  // dxdotdw[56]
    ddx_4acdt_dflux_re29 = 1;  // dxdotdw[57]
    ddx_k5k8k16dt_dflux_re30 = -1;  // dxdotdw[58]
    ddx_4acdt_dflux_re30 = 1;  // dxdotdw[59]
    ddx_k5k12k16dt_dflux_re31 = -1;  // dxdotdw[60]
    ddx_4acdt_dflux_re31 = 1;  // dxdotdw[61]
    ddx_k8k12k16dt_dflux_re32 = -1;  // dxdotdw[62]
    ddx_4acdt_dflux_re32 = 1;  // dxdotdw[63]
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
