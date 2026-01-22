#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 25> dwdp_colptrs_Blasi_CellSystems2016_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32
};

void dwdp_colptrs_Blasi_CellSystems2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Blasi_CellSystems2016_));
}
} // namespace model_Blasi_CellSystems2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<sunindextype, 32> dwdp_rowvals_Blasi_CellSystems2016_ = {
    26, 30, 35, 40, 43, 52, 56, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24
};

void dwdp_rowvals_Blasi_CellSystems2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Blasi_CellSystems2016_));
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
#include "dwdp.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void dwdp_Blasi_CellSystems2016(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        da_k12_da_basal = 1;  // dwdp[7]
        da_k12_k12k16_da_basal = 1;  // dwdp[8]
        da_k12_k8k12_da_basal = 1;  // dwdp[9]
        da_k12k16_k5k12k16_da_basal = 1;  // dwdp[10]
        da_k16_da_basal = 1;  // dwdp[11]
        da_k16_k5k16_da_basal = 1;  // dwdp[12]
        da_k16_k8k16_da_basal = 1;  // dwdp[13]
        da_k5_da_basal = 1;  // dwdp[14]
        da_k5_k5k16_da_basal = 1;  // dwdp[15]
        da_k5_k5k8_da_basal = 1;  // dwdp[16]
        da_k5k12_k5k12k16_da_basal = 1;  // dwdp[17]
        da_k5k12k16_4ac_da_basal = 1;  // dwdp[18]
        da_k5k16_k5k12k16_da_basal = 1;  // dwdp[19]
        da_k5k16_k5k8k16_da_basal = 1;  // dwdp[20]
        da_k5k8_k5k8k12_da_basal = 1;  // dwdp[21]
        da_k5k8_k5k8k16_da_basal = 1;  // dwdp[22]
        da_k5k8k12_4ac_da_basal = 1;  // dwdp[23]
        da_k5k8k16_4ac_da_basal = 1;  // dwdp[24]
        da_k8_k5k8_da_basal = 1;  // dwdp[25]
        da_k8_k8k12_da_basal = 1;  // dwdp[26]
        da_k8_k8k16_da_basal = 1;  // dwdp[27]
        da_k8k12_k5k8k12_da_basal = 1;  // dwdp[28]
        da_k8k12_k8k12k16_da_basal = 1;  // dwdp[29]
        da_k8k16_k5k8k16_da_basal = 1;  // dwdp[30]
        da_k8k16_k8k12k16_da_basal = 1;  // dwdp[31]
    }

    // dynamic expressions
    dflux_re2_da_k8 = x_0ac;  // dwdp[0]
    dflux_re6_da_k5_k5k12 = x_k5;  // dwdp[1]
    dflux_re11_da_k12_k5k12 = x_k12;  // dwdp[2]
    dflux_re16_da_k16_k12k16 = x_k16;  // dwdp[3]
    dflux_re19_da_k5k12_k5k8k12 = x_k5k12;  // dwdp[4]
    dflux_re28_da_k12k16_k8k12k16 = x_k12k16;  // dwdp[5]
    dflux_re32_da_k8k12k16_4ac = x_k8k12k16;  // dwdp[6]
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
