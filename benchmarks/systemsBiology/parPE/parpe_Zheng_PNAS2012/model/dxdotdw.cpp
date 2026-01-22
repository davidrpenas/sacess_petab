#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 62> dxdotdw_colptrs_Zheng_PNAS2012_ = {
    0, 0, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104
};

void dxdotdw_colptrs_Zheng_PNAS2012(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Zheng_PNAS2012_));
}
} // namespace model_Zheng_PNAS2012
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 104> dxdotdw_rowvals_Zheng_PNAS2012_ = {
    0, 0, 1, 0, 1, 0, 2, 0, 2, 2, 4, 2, 4, 1, 4, 1, 4, 1, 3, 1, 3, 2, 5, 2, 5, 3, 7, 3, 7, 5, 8, 5, 8, 3, 6, 3, 6, 5, 9, 5, 9, 4, 7, 4, 7, 4, 8, 4, 8, 7, 11, 7, 11, 8, 11, 8, 11, 11, 13, 11, 13, 11, 14, 11, 14, 9, 12, 9, 12, 6, 10, 6, 10, 7, 10, 7, 10, 8, 12, 8, 12, 12, 14, 12, 14, 10, 13, 10, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14
};

void dxdotdw_rowvals_Zheng_PNAS2012(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Zheng_PNAS2012_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void dxdotdw_Zheng_PNAS2012(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddK27me0K36me0dt_dflux_v1_v_0 = 1;  // dxdotdw[0]
    ddK27me0K36me0dt_dflux_v2_v_1 = -1;  // dxdotdw[1]
    ddK27me0K36me1dt_dflux_v2_v_1 = 1;  // dxdotdw[2]
    ddK27me0K36me0dt_dflux_v3_v_2 = 1;  // dxdotdw[3]
    ddK27me0K36me1dt_dflux_v3_v_2 = -1;  // dxdotdw[4]
    ddK27me0K36me0dt_dflux_v4_v_3 = -1;  // dxdotdw[5]
    ddK27me1K36me0dt_dflux_v4_v_3 = 1;  // dxdotdw[6]
    ddK27me0K36me0dt_dflux_v5_v_4 = 1;  // dxdotdw[7]
    ddK27me1K36me0dt_dflux_v5_v_4 = -1;  // dxdotdw[8]
    ddK27me1K36me0dt_dflux_v6_v_5 = -1;  // dxdotdw[9]
    ddK27me1K36me1dt_dflux_v6_v_5 = 1;  // dxdotdw[10]
    ddK27me1K36me0dt_dflux_v7_v_6 = 1;  // dxdotdw[11]
    ddK27me1K36me1dt_dflux_v7_v_6 = -1;  // dxdotdw[12]
    ddK27me0K36me1dt_dflux_v8_v_7 = -1;  // dxdotdw[13]
    ddK27me1K36me1dt_dflux_v8_v_7 = 1;  // dxdotdw[14]
    ddK27me0K36me1dt_dflux_v9_v_8 = 1;  // dxdotdw[15]
    ddK27me1K36me1dt_dflux_v9_v_8 = -1;  // dxdotdw[16]
    ddK27me0K36me1dt_dflux_v10_v_9 = -1;  // dxdotdw[17]
    ddK27me0K36me2dt_dflux_v10_v_9 = 1;  // dxdotdw[18]
    ddK27me0K36me1dt_dflux_v11_v_10 = 1;  // dxdotdw[19]
    ddK27me0K36me2dt_dflux_v11_v_10 = -1;  // dxdotdw[20]
    ddK27me1K36me0dt_dflux_v12_v_11 = -1;  // dxdotdw[21]
    ddK27me2K36me0dt_dflux_v12_v_11 = 1;  // dxdotdw[22]
    ddK27me1K36me0dt_dflux_v13_v_12 = 1;  // dxdotdw[23]
    ddK27me2K36me0dt_dflux_v13_v_12 = -1;  // dxdotdw[24]
    ddK27me0K36me2dt_dflux_v14_v_13 = -1;  // dxdotdw[25]
    ddK27me1K36me2dt_dflux_v14_v_13 = 1;  // dxdotdw[26]
    ddK27me0K36me2dt_dflux_v15_v_14 = 1;  // dxdotdw[27]
    ddK27me1K36me2dt_dflux_v15_v_14 = -1;  // dxdotdw[28]
    ddK27me2K36me0dt_dflux_v16_v_15 = -1;  // dxdotdw[29]
    ddK27me2K36me1dt_dflux_v16_v_15 = 1;  // dxdotdw[30]
    ddK27me2K36me0dt_dflux_v17_v_16 = 1;  // dxdotdw[31]
    ddK27me2K36me1dt_dflux_v17_v_16 = -1;  // dxdotdw[32]
    ddK27me0K36me2dt_dflux_v18_v_17 = -1;  // dxdotdw[33]
    ddK27me0K36me3dt_dflux_v18_v_17 = 1;  // dxdotdw[34]
    ddK27me0K36me2dt_dflux_v19_v_18 = 1;  // dxdotdw[35]
    ddK27me0K36me3dt_dflux_v19_v_18 = -1;  // dxdotdw[36]
    ddK27me2K36me0dt_dflux_v20_v_19 = -1;  // dxdotdw[37]
    ddK27me3K36me0dt_dflux_v20_v_19 = 1;  // dxdotdw[38]
    ddK27me2K36me0dt_dflux_v21_v_20 = 1;  // dxdotdw[39]
    ddK27me3K36me0dt_dflux_v21_v_20 = -1;  // dxdotdw[40]
    ddK27me1K36me1dt_dflux_v22_v_21 = -1;  // dxdotdw[41]
    ddK27me1K36me2dt_dflux_v22_v_21 = 1;  // dxdotdw[42]
    ddK27me1K36me1dt_dflux_v23_v_22 = 1;  // dxdotdw[43]
    ddK27me1K36me2dt_dflux_v23_v_22 = -1;  // dxdotdw[44]
    ddK27me1K36me1dt_dflux_v24_v_23 = -1;  // dxdotdw[45]
    ddK27me2K36me1dt_dflux_v24_v_23 = 1;  // dxdotdw[46]
    ddK27me1K36me1dt_dflux_v25_v_24 = 1;  // dxdotdw[47]
    ddK27me2K36me1dt_dflux_v25_v_24 = -1;  // dxdotdw[48]
    ddK27me1K36me2dt_dflux_v26_v_25 = -1;  // dxdotdw[49]
    ddK27me2K36me2dt_dflux_v26_v_25 = 1;  // dxdotdw[50]
    ddK27me1K36me2dt_dflux_v27_v_26 = 1;  // dxdotdw[51]
    ddK27me2K36me2dt_dflux_v27_v_26 = -1;  // dxdotdw[52]
    ddK27me2K36me1dt_dflux_v28_v_27 = -1;  // dxdotdw[53]
    ddK27me2K36me2dt_dflux_v28_v_27 = 1;  // dxdotdw[54]
    ddK27me2K36me1dt_dflux_v29_v_28 = 1;  // dxdotdw[55]
    ddK27me2K36me2dt_dflux_v29_v_28 = -1;  // dxdotdw[56]
    ddK27me2K36me2dt_dflux_v30_v_29 = -1;  // dxdotdw[57]
    ddK27me2K36me3dt_dflux_v30_v_29 = 1;  // dxdotdw[58]
    ddK27me2K36me2dt_dflux_v31_v_30 = 1;  // dxdotdw[59]
    ddK27me2K36me3dt_dflux_v31_v_30 = -1;  // dxdotdw[60]
    ddK27me2K36me2dt_dflux_v32_v_31 = -1;  // dxdotdw[61]
    ddK27me3K36me2dt_dflux_v32_v_31 = 1;  // dxdotdw[62]
    ddK27me2K36me2dt_dflux_v33_v_32 = 1;  // dxdotdw[63]
    ddK27me3K36me2dt_dflux_v33_v_32 = -1;  // dxdotdw[64]
    ddK27me3K36me0dt_dflux_v34_v_33 = -1;  // dxdotdw[65]
    ddK27me3K36me1dt_dflux_v34_v_33 = 1;  // dxdotdw[66]
    ddK27me3K36me0dt_dflux_v35_v_34 = 1;  // dxdotdw[67]
    ddK27me3K36me1dt_dflux_v35_v_34 = -1;  // dxdotdw[68]
    ddK27me0K36me3dt_dflux_v36_v_35 = -1;  // dxdotdw[69]
    ddK27me1K36me3dt_dflux_v36_v_35 = 1;  // dxdotdw[70]
    ddK27me0K36me3dt_dflux_v37_v_36 = 1;  // dxdotdw[71]
    ddK27me1K36me3dt_dflux_v37_v_36 = -1;  // dxdotdw[72]
    ddK27me1K36me2dt_dflux_v38_v_37 = -1;  // dxdotdw[73]
    ddK27me1K36me3dt_dflux_v38_v_37 = 1;  // dxdotdw[74]
    ddK27me1K36me2dt_dflux_v39_v_38 = 1;  // dxdotdw[75]
    ddK27me1K36me3dt_dflux_v39_v_38 = -1;  // dxdotdw[76]
    ddK27me2K36me1dt_dflux_v40_v_39 = -1;  // dxdotdw[77]
    ddK27me3K36me1dt_dflux_v40_v_39 = 1;  // dxdotdw[78]
    ddK27me2K36me1dt_dflux_v41_v_40 = 1;  // dxdotdw[79]
    ddK27me3K36me1dt_dflux_v41_v_40 = -1;  // dxdotdw[80]
    ddK27me3K36me1dt_dflux_v42_v_41 = -1;  // dxdotdw[81]
    ddK27me3K36me2dt_dflux_v42_v_41 = 1;  // dxdotdw[82]
    ddK27me3K36me1dt_dflux_v43_v_42 = 1;  // dxdotdw[83]
    ddK27me3K36me2dt_dflux_v43_v_42 = -1;  // dxdotdw[84]
    ddK27me1K36me3dt_dflux_v44_v_43 = -1;  // dxdotdw[85]
    ddK27me2K36me3dt_dflux_v44_v_43 = 1;  // dxdotdw[86]
    ddK27me1K36me3dt_dflux_v45_v_44 = 1;  // dxdotdw[87]
    ddK27me2K36me3dt_dflux_v45_v_44 = -1;  // dxdotdw[88]
    ddK27me0K36me0dt_dflux_v46_v_45 = -1;  // dxdotdw[89]
    ddK27me0K36me1dt_dflux_v47_v_46 = -1;  // dxdotdw[90]
    ddK27me1K36me0dt_dflux_v48_v_47 = -1;  // dxdotdw[91]
    ddK27me0K36me2dt_dflux_v49_v_48 = -1;  // dxdotdw[92]
    ddK27me1K36me1dt_dflux_v50_v_49 = -1;  // dxdotdw[93]
    ddK27me2K36me0dt_dflux_v51_v_50 = -1;  // dxdotdw[94]
    ddK27me0K36me3dt_dflux_v52_v_51 = -1;  // dxdotdw[95]
    ddK27me1K36me2dt_dflux_v53_v_52 = -1;  // dxdotdw[96]
    ddK27me2K36me1dt_dflux_v54_v_53 = -1;  // dxdotdw[97]
    ddK27me3K36me0dt_dflux_v55_v_54 = -1;  // dxdotdw[98]
    ddK27me1K36me3dt_dflux_v56_v_55 = -1;  // dxdotdw[99]
    ddK27me2K36me2dt_dflux_v57_v_56 = -1;  // dxdotdw[100]
    ddK27me3K36me1dt_dflux_v58_v_57 = -1;  // dxdotdw[101]
    ddK27me2K36me3dt_dflux_v59_v_58 = -1;  // dxdotdw[102]
    ddK27me3K36me2dt_dflux_v60_v_59 = -1;  // dxdotdw[103]
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
