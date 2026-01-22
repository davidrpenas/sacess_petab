#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 17> dwdx_colptrs_Blasi_CellSystems2016_ = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64
};

void dwdx_colptrs_Blasi_CellSystems2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Blasi_CellSystems2016_));
}
} // namespace model_Blasi_CellSystems2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 64> dwdx_rowvals_Blasi_CellSystems2016_ = {
    25, 26, 27, 28, 26, 32, 33, 34, 25, 29, 30, 31, 27, 35, 36, 37, 28, 38, 39, 40, 29, 32, 41, 42, 30, 35, 43, 44, 31, 38, 45, 46, 33, 36, 47, 48, 34, 39, 49, 50, 37, 40, 51, 52, 41, 43, 47, 53, 42, 45, 49, 54, 44, 46, 51, 55, 48, 50, 52, 56, 53, 54, 55, 56
};

void dwdx_rowvals_Blasi_CellSystems2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Blasi_CellSystems2016_));
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
#include "dwdx.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void dwdx_Blasi_CellSystems2016(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_re1_dx_0ac = a_k5;  // dwdx[0]
        dflux_re2_dx_0ac = a_k8;  // dwdx[1]
        dflux_re3_dx_0ac = a_k12;  // dwdx[2]
        dflux_re4_dx_0ac = a_k16;  // dwdx[3]
        dflux_re2_dx_k8 = -d;  // dwdx[4]
        dflux_re8_dx_k8 = a_k8_k5k8;  // dwdx[5]
        dflux_re9_dx_k8 = a_k8_k8k12;  // dwdx[6]
        dflux_re10_dx_k8 = a_k8_k8k16;  // dwdx[7]
        dflux_re1_dx_k5 = -d;  // dwdx[8]
        dflux_re5_dx_k5 = a_k5_k5k8;  // dwdx[9]
        dflux_re6_dx_k5 = a_k5_k5k12;  // dwdx[10]
        dflux_re7_dx_k5 = a_k5_k5k16;  // dwdx[11]
        dflux_re3_dx_k12 = -d;  // dwdx[12]
        dflux_re11_dx_k12 = a_k12_k5k12;  // dwdx[13]
        dflux_re12_dx_k12 = a_k12_k8k12;  // dwdx[14]
        dflux_re13_dx_k12 = a_k12_k12k16;  // dwdx[15]
        dflux_re4_dx_k16 = -d;  // dwdx[16]
        dflux_re14_dx_k16 = a_k16_k5k16;  // dwdx[17]
        dflux_re15_dx_k16 = a_k16_k8k16;  // dwdx[18]
        dflux_re16_dx_k16 = a_k16_k12k16;  // dwdx[19]
        dflux_re5_dx_k5k8 = -d;  // dwdx[20]
        dflux_re8_dx_k5k8 = -d;  // dwdx[21]
        dflux_re17_dx_k5k8 = a_k5k8_k5k8k12;  // dwdx[22]
        dflux_re18_dx_k5k8 = a_k5k8_k5k8k16;  // dwdx[23]
        dflux_re6_dx_k5k12 = -d;  // dwdx[24]
        dflux_re11_dx_k5k12 = -d;  // dwdx[25]
        dflux_re19_dx_k5k12 = a_k5k12_k5k8k12;  // dwdx[26]
        dflux_re20_dx_k5k12 = a_k5k12_k5k12k16;  // dwdx[27]
        dflux_re7_dx_k5k16 = -d;  // dwdx[28]
        dflux_re14_dx_k5k16 = -d;  // dwdx[29]
        dflux_re21_dx_k5k16 = a_k5k16_k5k8k16;  // dwdx[30]
        dflux_re22_dx_k5k16 = a_k5k16_k5k12k16;  // dwdx[31]
        dflux_re9_dx_k8k12 = -d;  // dwdx[32]
        dflux_re12_dx_k8k12 = -d;  // dwdx[33]
        dflux_re23_dx_k8k12 = a_k8k12_k5k8k12;  // dwdx[34]
        dflux_re24_dx_k8k12 = a_k8k12_k8k12k16;  // dwdx[35]
        dflux_re10_dx_k8k16 = -d;  // dwdx[36]
        dflux_re15_dx_k8k16 = -d;  // dwdx[37]
        dflux_re25_dx_k8k16 = a_k8k16_k5k8k16;  // dwdx[38]
        dflux_re26_dx_k8k16 = a_k8k16_k8k12k16;  // dwdx[39]
        dflux_re13_dx_k12k16 = -d;  // dwdx[40]
        dflux_re16_dx_k12k16 = -d;  // dwdx[41]
        dflux_re27_dx_k12k16 = a_k12k16_k5k12k16;  // dwdx[42]
        dflux_re28_dx_k12k16 = a_k12k16_k8k12k16;  // dwdx[43]
        dflux_re17_dx_k5k8k12 = -d;  // dwdx[44]
        dflux_re19_dx_k5k8k12 = -d;  // dwdx[45]
        dflux_re23_dx_k5k8k12 = -d;  // dwdx[46]
        dflux_re29_dx_k5k8k12 = a_k5k8k12_4ac;  // dwdx[47]
        dflux_re18_dx_k5k8k16 = -d;  // dwdx[48]
        dflux_re21_dx_k5k8k16 = -d;  // dwdx[49]
        dflux_re25_dx_k5k8k16 = -d;  // dwdx[50]
        dflux_re30_dx_k5k8k16 = a_k5k8k16_4ac;  // dwdx[51]
        dflux_re20_dx_k5k12k16 = -d;  // dwdx[52]
        dflux_re22_dx_k5k12k16 = -d;  // dwdx[53]
        dflux_re27_dx_k5k12k16 = -d;  // dwdx[54]
        dflux_re31_dx_k5k12k16 = a_k5k12k16_4ac;  // dwdx[55]
        dflux_re24_dx_k8k12k16 = -d;  // dwdx[56]
        dflux_re26_dx_k8k12k16 = -d;  // dwdx[57]
        dflux_re28_dx_k8k12k16 = -d;  // dwdx[58]
        dflux_re32_dx_k8k12k16 = a_k8k12k16_4ac;  // dwdx[59]
        dflux_re29_dx_4ac = -d;  // dwdx[60]
        dflux_re30_dx_4ac = -d;  // dwdx[61]
        dflux_re31_dx_4ac = -d;  // dwdx[62]
        dflux_re32_dx_4ac = -d;  // dwdx[63]
    }
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
