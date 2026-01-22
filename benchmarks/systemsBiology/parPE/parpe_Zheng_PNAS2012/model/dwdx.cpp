#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 16> dwdx_colptrs_Zheng_PNAS2012_ = {
    0, 3, 7, 11, 15, 20, 24, 27, 32, 37, 40, 44, 49, 53, 56, 59
};

void dwdx_colptrs_Zheng_PNAS2012(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Zheng_PNAS2012_));
}
} // namespace model_Zheng_PNAS2012
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 59> dwdx_rowvals_Zheng_PNAS2012_ = {
    2, 4, 46, 3, 8, 10, 47, 5, 6, 12, 48, 11, 14, 18, 49, 7, 9, 22, 24, 50, 13, 16, 20, 51, 19, 36, 52, 15, 23, 26, 38, 53, 17, 25, 28, 40, 54, 21, 34, 55, 37, 39, 44, 56, 27, 29, 30, 32, 57, 35, 41, 42, 58, 31, 45, 59, 33, 43, 60
};

void dwdx_rowvals_Zheng_PNAS2012(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Zheng_PNAS2012_));
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
#include "dwdx.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void dwdx_Zheng_PNAS2012(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v2_v_1_dK27me0K36me0 = k00_01;  // dwdx[0]
        dflux_v4_v_3_dK27me0K36me0 = k00_10;  // dwdx[1]
        dflux_v46_v_45_dK27me0K36me0 = inflow;  // dwdx[2]
        dflux_v3_v_2_dK27me0K36me1 = k01_00;  // dwdx[3]
        dflux_v8_v_7_dK27me0K36me1 = k01_11;  // dwdx[4]
        dflux_v10_v_9_dK27me0K36me1 = k01_02;  // dwdx[5]
        dflux_v47_v_46_dK27me0K36me1 = inflow;  // dwdx[6]
        dflux_v5_v_4_dK27me1K36me0 = k10_00;  // dwdx[7]
        dflux_v6_v_5_dK27me1K36me0 = k10_11;  // dwdx[8]
        dflux_v12_v_11_dK27me1K36me0 = k10_20;  // dwdx[9]
        dflux_v48_v_47_dK27me1K36me0 = inflow;  // dwdx[10]
        dflux_v11_v_10_dK27me0K36me2 = k02_01;  // dwdx[11]
        dflux_v14_v_13_dK27me0K36me2 = k02_12;  // dwdx[12]
        dflux_v18_v_17_dK27me0K36me2 = k02_03;  // dwdx[13]
        dflux_v49_v_48_dK27me0K36me2 = inflow;  // dwdx[14]
        dflux_v7_v_6_dK27me1K36me1 = k11_10;  // dwdx[15]
        dflux_v9_v_8_dK27me1K36me1 = k11_01;  // dwdx[16]
        dflux_v22_v_21_dK27me1K36me1 = k11_12;  // dwdx[17]
        dflux_v24_v_23_dK27me1K36me1 = k11_21;  // dwdx[18]
        dflux_v50_v_49_dK27me1K36me1 = inflow;  // dwdx[19]
        dflux_v13_v_12_dK27me2K36me0 = k20_10;  // dwdx[20]
        dflux_v16_v_15_dK27me2K36me0 = k20_21;  // dwdx[21]
        dflux_v20_v_19_dK27me2K36me0 = k20_30;  // dwdx[22]
        dflux_v51_v_50_dK27me2K36me0 = inflow;  // dwdx[23]
        dflux_v19_v_18_dK27me0K36me3 = k03_02;  // dwdx[24]
        dflux_v36_v_35_dK27me0K36me3 = k03_13;  // dwdx[25]
        dflux_v52_v_51_dK27me0K36me3 = inflow;  // dwdx[26]
        dflux_v15_v_14_dK27me1K36me2 = k12_02;  // dwdx[27]
        dflux_v23_v_22_dK27me1K36me2 = k12_11;  // dwdx[28]
        dflux_v26_v_25_dK27me1K36me2 = k12_22;  // dwdx[29]
        dflux_v38_v_37_dK27me1K36me2 = k12_13;  // dwdx[30]
        dflux_v53_v_52_dK27me1K36me2 = inflow;  // dwdx[31]
        dflux_v17_v_16_dK27me2K36me1 = k21_20;  // dwdx[32]
        dflux_v25_v_24_dK27me2K36me1 = k21_11;  // dwdx[33]
        dflux_v28_v_27_dK27me2K36me1 = k21_22;  // dwdx[34]
        dflux_v40_v_39_dK27me2K36me1 = k21_31;  // dwdx[35]
        dflux_v54_v_53_dK27me2K36me1 = inflow;  // dwdx[36]
        dflux_v21_v_20_dK27me3K36me0 = k30_20;  // dwdx[37]
        dflux_v34_v_33_dK27me3K36me0 = k30_31;  // dwdx[38]
        dflux_v55_v_54_dK27me3K36me0 = inflow;  // dwdx[39]
        dflux_v37_v_36_dK27me1K36me3 = k13_03;  // dwdx[40]
        dflux_v39_v_38_dK27me1K36me3 = k13_12;  // dwdx[41]
        dflux_v44_v_43_dK27me1K36me3 = k13_23;  // dwdx[42]
        dflux_v56_v_55_dK27me1K36me3 = inflow;  // dwdx[43]
        dflux_v27_v_26_dK27me2K36me2 = k22_12;  // dwdx[44]
        dflux_v29_v_28_dK27me2K36me2 = k22_21;  // dwdx[45]
        dflux_v30_v_29_dK27me2K36me2 = k22_23;  // dwdx[46]
        dflux_v32_v_31_dK27me2K36me2 = k22_32;  // dwdx[47]
        dflux_v57_v_56_dK27me2K36me2 = inflow;  // dwdx[48]
        dflux_v35_v_34_dK27me3K36me1 = k31_30;  // dwdx[49]
        dflux_v41_v_40_dK27me3K36me1 = k31_21;  // dwdx[50]
        dflux_v42_v_41_dK27me3K36me1 = k31_32;  // dwdx[51]
        dflux_v58_v_57_dK27me3K36me1 = inflow;  // dwdx[52]
        dflux_v31_v_30_dK27me2K36me3 = k23_22;  // dwdx[53]
        dflux_v45_v_44_dK27me2K36me3 = k23_13;  // dwdx[54]
        dflux_v59_v_58_dK27me2K36me3 = inflow;  // dwdx[55]
        dflux_v33_v_32_dK27me3K36me2 = k32_22;  // dwdx[56]
        dflux_v43_v_42_dK27me3K36me2 = k32_31;  // dwdx[57]
        dflux_v60_v_59_dK27me3K36me2 = inflow;  // dwdx[58]
    }
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
