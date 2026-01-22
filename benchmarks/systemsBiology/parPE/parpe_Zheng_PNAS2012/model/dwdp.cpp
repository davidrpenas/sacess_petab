#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 61> dwdp_colptrs_Zheng_PNAS2012_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45
};

void dwdp_colptrs_Zheng_PNAS2012(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Zheng_PNAS2012_));
}
} // namespace model_Zheng_PNAS2012
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 45> dwdp_rowvals_Zheng_PNAS2012_ = {
    0, 2, 4, 3, 10, 8, 11, 18, 14, 19, 36, 5, 6, 12, 9, 7, 22, 24, 15, 23, 38, 26, 37, 39, 44, 13, 16, 20, 25, 17, 28, 40, 27, 29, 30, 32, 45, 31, 21, 34, 41, 35, 42, 33, 43
};

void dwdp_rowvals_Zheng_PNAS2012(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Zheng_PNAS2012_));
}
} // namespace model_Zheng_PNAS2012
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
namespace model_Zheng_PNAS2012 {

void dwdp_Zheng_PNAS2012(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dinflow_dinflowp = dilution;  // dwdp[0]
    }

    // dynamic expressions
    dflux_v2_v_1_dk00_01 = K27me0K36me0;  // dwdp[1]
    dflux_v4_v_3_dk00_10 = K27me0K36me0;  // dwdp[2]
    dflux_v3_v_2_dk01_00 = K27me0K36me1;  // dwdp[3]
    dflux_v10_v_9_dk01_02 = K27me0K36me1;  // dwdp[4]
    dflux_v8_v_7_dk01_11 = K27me0K36me1;  // dwdp[5]
    dflux_v11_v_10_dk02_01 = K27me0K36me2;  // dwdp[6]
    dflux_v18_v_17_dk02_03 = K27me0K36me2;  // dwdp[7]
    dflux_v14_v_13_dk02_12 = K27me0K36me2;  // dwdp[8]
    dflux_v19_v_18_dk03_02 = K27me0K36me3;  // dwdp[9]
    dflux_v36_v_35_dk03_13 = K27me0K36me3;  // dwdp[10]
    dflux_v5_v_4_dk10_00 = K27me1K36me0;  // dwdp[11]
    dflux_v6_v_5_dk10_11 = K27me1K36me0;  // dwdp[12]
    dflux_v12_v_11_dk10_20 = K27me1K36me0;  // dwdp[13]
    dflux_v9_v_8_dk11_01 = K27me1K36me1;  // dwdp[14]
    dflux_v7_v_6_dk11_10 = K27me1K36me1;  // dwdp[15]
    dflux_v22_v_21_dk11_12 = K27me1K36me1;  // dwdp[16]
    dflux_v24_v_23_dk11_21 = K27me1K36me1;  // dwdp[17]
    dflux_v15_v_14_dk12_02 = K27me1K36me2;  // dwdp[18]
    dflux_v23_v_22_dk12_11 = K27me1K36me2;  // dwdp[19]
    dflux_v38_v_37_dk12_13 = K27me1K36me2;  // dwdp[20]
    dflux_v26_v_25_dk12_22 = K27me1K36me2;  // dwdp[21]
    dflux_v37_v_36_dk13_03 = K27me1K36me3;  // dwdp[22]
    dflux_v39_v_38_dk13_12 = K27me1K36me3;  // dwdp[23]
    dflux_v44_v_43_dk13_23 = K27me1K36me3;  // dwdp[24]
    dflux_v13_v_12_dk20_10 = K27me2K36me0;  // dwdp[25]
    dflux_v16_v_15_dk20_21 = K27me2K36me0;  // dwdp[26]
    dflux_v20_v_19_dk20_30 = K27me2K36me0;  // dwdp[27]
    dflux_v25_v_24_dk21_11 = K27me2K36me1;  // dwdp[28]
    dflux_v17_v_16_dk21_20 = K27me2K36me1;  // dwdp[29]
    dflux_v28_v_27_dk21_22 = K27me2K36me1;  // dwdp[30]
    dflux_v40_v_39_dk21_31 = K27me2K36me1;  // dwdp[31]
    dflux_v27_v_26_dk22_12 = K27me2K36me2;  // dwdp[32]
    dflux_v29_v_28_dk22_21 = K27me2K36me2;  // dwdp[33]
    dflux_v30_v_29_dk22_23 = K27me2K36me2;  // dwdp[34]
    dflux_v32_v_31_dk22_32 = K27me2K36me2;  // dwdp[35]
    dflux_v45_v_44_dk23_13 = K27me2K36me3;  // dwdp[36]
    dflux_v31_v_30_dk23_22 = K27me2K36me3;  // dwdp[37]
    dflux_v21_v_20_dk30_20 = K27me3K36me0;  // dwdp[38]
    dflux_v34_v_33_dk30_31 = K27me3K36me0;  // dwdp[39]
    dflux_v41_v_40_dk31_21 = K27me3K36me1;  // dwdp[40]
    dflux_v35_v_34_dk31_30 = K27me3K36me1;  // dwdp[41]
    dflux_v42_v_41_dk31_32 = K27me3K36me1;  // dwdp[42]
    dflux_v33_v_32_dk32_22 = K27me3K36me2;  // dwdp[43]
    dflux_v43_v_42_dk32_31 = K27me3K36me2;  // dwdp[44]
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
