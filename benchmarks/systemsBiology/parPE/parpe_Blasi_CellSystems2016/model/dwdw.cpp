#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 58> dwdw_colptrs_Blasi_CellSystems2016_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25
};

void dwdw_colptrs_Blasi_CellSystems2016(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Blasi_CellSystems2016_));
}
} // namespace model_Blasi_CellSystems2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 25> dwdw_rowvals_Blasi_CellSystems2016_ = {
    27, 37, 36, 51, 28, 38, 39, 25, 31, 29, 44, 55, 46, 45, 41, 42, 53, 54, 32, 33, 34, 47, 48, 49, 50
};

void dwdw_rowvals_Blasi_CellSystems2016(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Blasi_CellSystems2016_));
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
#include "dwdw.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void dwdw_Blasi_CellSystems2016(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_re3_da_k12 = x_0ac;  // dwdw[0]
    dflux_re13_da_k12_k12k16 = x_k12;  // dwdw[1]
    dflux_re12_da_k12_k8k12 = x_k12;  // dwdw[2]
    dflux_re27_da_k12k16_k5k12k16 = x_k12k16;  // dwdw[3]
    dflux_re4_da_k16 = x_0ac;  // dwdw[4]
    dflux_re14_da_k16_k5k16 = x_k16;  // dwdw[5]
    dflux_re15_da_k16_k8k16 = x_k16;  // dwdw[6]
    dflux_re1_da_k5 = x_0ac;  // dwdw[7]
    dflux_re7_da_k5_k5k16 = x_k5;  // dwdw[8]
    dflux_re5_da_k5_k5k8 = x_k5;  // dwdw[9]
    dflux_re20_da_k5k12_k5k12k16 = x_k5k12;  // dwdw[10]
    dflux_re31_da_k5k12k16_4ac = x_k5k12k16;  // dwdw[11]
    dflux_re22_da_k5k16_k5k12k16 = x_k5k16;  // dwdw[12]
    dflux_re21_da_k5k16_k5k8k16 = x_k5k16;  // dwdw[13]
    dflux_re17_da_k5k8_k5k8k12 = x_k5k8;  // dwdw[14]
    dflux_re18_da_k5k8_k5k8k16 = x_k5k8;  // dwdw[15]
    dflux_re29_da_k5k8k12_4ac = x_k5k8k12;  // dwdw[16]
    dflux_re30_da_k5k8k16_4ac = x_k5k8k16;  // dwdw[17]
    dflux_re8_da_k8_k5k8 = x_k8;  // dwdw[18]
    dflux_re9_da_k8_k8k12 = x_k8;  // dwdw[19]
    dflux_re10_da_k8_k8k16 = x_k8;  // dwdw[20]
    dflux_re23_da_k8k12_k5k8k12 = x_k8k12;  // dwdw[21]
    dflux_re24_da_k8k12_k8k12k16 = x_k8k12;  // dwdw[22]
    dflux_re25_da_k8k16_k5k8k16 = x_k8k16;  // dwdw[23]
    dflux_re26_da_k8k16_k8k12k16 = x_k8k16;  // dwdw[24]
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
