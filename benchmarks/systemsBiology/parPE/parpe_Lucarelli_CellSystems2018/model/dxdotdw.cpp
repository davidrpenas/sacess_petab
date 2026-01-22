#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<sunindextype, 58> dxdotdw_colptrs_Lucarelli_CellSystems2018_ = {
    0, 3, 4, 7, 10, 12, 15, 17, 20, 22, 24, 26, 28, 30, 32, 34, 36, 39, 43, 46, 49, 53, 56, 59, 63, 66, 70, 73, 76, 79, 82, 86, 90, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118
};

void dxdotdw_colptrs_Lucarelli_CellSystems2018(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Lucarelli_CellSystems2018_));
}
} // namespace model_Lucarelli_CellSystems2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

static constexpr std::array<sunindextype, 118> dxdotdw_rowvals_Lucarelli_CellSystems2018_ = {
    0, 1, 2, 2, 3, 5, 6, 3, 5, 6, 7, 12, 7, 10, 12, 8, 13, 8, 11, 13, 5, 9, 5, 9, 3, 12, 10, 12, 3, 10, 4, 13, 11, 13, 4, 11, 12, 13, 15, 10, 12, 13, 15, 11, 12, 15, 5, 12, 14, 5, 10, 12, 14, 12, 13, 16, 10, 13, 16, 11, 12, 13, 16, 5, 13, 17, 5, 11, 13, 17, 5, 12, 20, 5, 10, 20, 5, 13, 19, 5, 11, 19, 5, 12, 13, 18, 5, 10, 13, 18, 5, 11, 12, 18, 21, 21, 22, 22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32
};

void dxdotdw_rowvals_Lucarelli_CellSystems2018(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Lucarelli_CellSystems2018_));
}
} // namespace model_Lucarelli_CellSystems2018
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
namespace model_Lucarelli_CellSystems2018 {

void dxdotdw_Lucarelli_CellSystems2018(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddTGFbdt_dflux_v1_v_0 = -0.075930144267274111;  // dxdotdw[0]
    ddRecdt_dflux_v1_v_0 = -0.075930144267274111;  // dxdotdw[1]
    ddTGFb_pRecdt_dflux_v1_v_0 = 0.075930144267274111;  // dxdotdw[2]
    ddTGFb_pRecdt_dflux_v2_v_1 = -0.075930144267274111;  // dxdotdw[3]
    ddS2dt_dflux_v3_v_2 = -0.075930144267274111;  // dxdotdw[4]
    ddS4dt_dflux_v3_v_2 = -0.15186028853454822;  // dxdotdw[5]
    ddS2_S4_S4dt_dflux_v3_v_2 = 0.075930144267274111;  // dxdotdw[6]
    ddS2dt_dflux_v4_v_3 = 0.075930144267274111;  // dxdotdw[7]
    ddS4dt_dflux_v4_v_3 = 0.15186028853454822;  // dxdotdw[8]
    ddS2_S4_S4dt_dflux_v4_v_3 = -0.075930144267274111;  // dxdotdw[9]
    ddppS2_ppS2_ppS2dt_dflux_v5_v_4 = 0.075930144267274111;  // dxdotdw[10]
    ddppS2dt_dflux_v5_v_4 = -0.22779043280182232;  // dxdotdw[11]
    ddppS2_ppS2_ppS2dt_dflux_v6_v_5 = -0.075930144267274111;  // dxdotdw[12]
    ddpS2dt_dflux_v6_v_5 = 0.075930144267274111;  // dxdotdw[13]
    ddppS2dt_dflux_v6_v_5 = 0.15186028853454822;  // dxdotdw[14]
    ddppS3_ppS3_ppS3dt_dflux_v7_v_6 = 0.075930144267274111;  // dxdotdw[15]
    ddppS3dt_dflux_v7_v_6 = -0.22779043280182232;  // dxdotdw[16]
    ddppS3_ppS3_ppS3dt_dflux_v8_v_7 = -0.075930144267274111;  // dxdotdw[17]
    ddpS3dt_dflux_v8_v_7 = 0.075930144267274111;  // dxdotdw[18]
    ddppS3dt_dflux_v8_v_7 = 0.15186028853454822;  // dxdotdw[19]
    ddS4dt_dflux_v9_v_8 = -0.22779043280182232;  // dxdotdw[20]
    ddS4_S4_S4dt_dflux_v9_v_8 = 0.075930144267274111;  // dxdotdw[21]
    ddS4dt_dflux_v10_v_9 = 0.22779043280182232;  // dxdotdw[22]
    ddS4_S4_S4dt_dflux_v10_v_9 = -0.075930144267274111;  // dxdotdw[23]
    ddS2dt_dflux_v11_v_10 = -0.075930144267274111;  // dxdotdw[24]
    ddppS2dt_dflux_v11_v_10 = 0.075930144267274111;  // dxdotdw[25]
    ddpS2dt_dflux_v12_v_11 = 0.075930144267274111;  // dxdotdw[26]
    ddppS2dt_dflux_v12_v_11 = -0.075930144267274111;  // dxdotdw[27]
    ddS2dt_dflux_v13_v_12 = 0.075930144267274111;  // dxdotdw[28]
    ddpS2dt_dflux_v13_v_12 = -0.075930144267274111;  // dxdotdw[29]
    ddS3dt_dflux_v14_v_13 = -0.075930144267274111;  // dxdotdw[30]
    ddppS3dt_dflux_v14_v_13 = 0.075930144267274111;  // dxdotdw[31]
    ddpS3dt_dflux_v15_v_14 = 0.075930144267274111;  // dxdotdw[32]
    ddppS3dt_dflux_v15_v_14 = -0.075930144267274111;  // dxdotdw[33]
    ddS3dt_dflux_v16_v_15 = 0.075930144267274111;  // dxdotdw[34]
    ddpS3dt_dflux_v16_v_15 = -0.075930144267274111;  // dxdotdw[35]
    ddppS2dt_dflux_v17_v_16 = -0.15186028853454822;  // dxdotdw[36]
    ddppS3dt_dflux_v17_v_16 = -0.075930144267274111;  // dxdotdw[37]
    ddppS2_ppS2_ppS3dt_dflux_v17_v_16 = 0.075930144267274111;  // dxdotdw[38]
    ddpS2dt_dflux_v18_v_17 = 0.075930144267274111;  // dxdotdw[39]
    ddppS2dt_dflux_v18_v_17 = 0.075930144267274111;  // dxdotdw[40]
    ddppS3dt_dflux_v18_v_17 = 0.075930144267274111;  // dxdotdw[41]
    ddppS2_ppS2_ppS3dt_dflux_v18_v_17 = -0.075930144267274111;  // dxdotdw[42]
    ddpS3dt_dflux_v19_v_18 = 0.075930144267274111;  // dxdotdw[43]
    ddppS2dt_dflux_v19_v_18 = 0.15186028853454822;  // dxdotdw[44]
    ddppS2_ppS2_ppS3dt_dflux_v19_v_18 = -0.075930144267274111;  // dxdotdw[45]
    ddS4dt_dflux_v20_v_19 = -0.075930144267274111;  // dxdotdw[46]
    ddppS2dt_dflux_v20_v_19 = -0.15186028853454822;  // dxdotdw[47]
    ddppS2_ppS2_S4dt_dflux_v20_v_19 = 0.075930144267274111;  // dxdotdw[48]
    ddS4dt_dflux_v21_v_20 = 0.075930144267274111;  // dxdotdw[49]
    ddpS2dt_dflux_v21_v_20 = 0.075930144267274111;  // dxdotdw[50]
    ddppS2dt_dflux_v21_v_20 = 0.075930144267274111;  // dxdotdw[51]
    ddppS2_ppS2_S4dt_dflux_v21_v_20 = -0.075930144267274111;  // dxdotdw[52]
    ddppS2dt_dflux_v22_v_21 = -0.075930144267274111;  // dxdotdw[53]
    ddppS3dt_dflux_v22_v_21 = -0.15186028853454822;  // dxdotdw[54]
    ddppS2_ppS3_ppS3dt_dflux_v22_v_21 = 0.075930144267274111;  // dxdotdw[55]
    ddpS2dt_dflux_v23_v_22 = 0.075930144267274111;  // dxdotdw[56]
    ddppS3dt_dflux_v23_v_22 = 0.15186028853454822;  // dxdotdw[57]
    ddppS2_ppS3_ppS3dt_dflux_v23_v_22 = -0.075930144267274111;  // dxdotdw[58]
    ddpS3dt_dflux_v24_v_23 = 0.075930144267274111;  // dxdotdw[59]
    ddppS2dt_dflux_v24_v_23 = 0.075930144267274111;  // dxdotdw[60]
    ddppS3dt_dflux_v24_v_23 = 0.075930144267274111;  // dxdotdw[61]
    ddppS2_ppS3_ppS3dt_dflux_v24_v_23 = -0.075930144267274111;  // dxdotdw[62]
    ddS4dt_dflux_v25_v_24 = -0.075930144267274111;  // dxdotdw[63]
    ddppS3dt_dflux_v25_v_24 = -0.15186028853454822;  // dxdotdw[64]
    ddppS3_ppS3_S4dt_dflux_v25_v_24 = 0.075930144267274111;  // dxdotdw[65]
    ddS4dt_dflux_v26_v_25 = 0.075930144267274111;  // dxdotdw[66]
    ddpS3dt_dflux_v26_v_25 = 0.075930144267274111;  // dxdotdw[67]
    ddppS3dt_dflux_v26_v_25 = 0.075930144267274111;  // dxdotdw[68]
    ddppS3_ppS3_S4dt_dflux_v26_v_25 = -0.075930144267274111;  // dxdotdw[69]
    ddS4dt_dflux_v27_v_26 = -0.15186028853454822;  // dxdotdw[70]
    ddppS2dt_dflux_v27_v_26 = -0.075930144267274111;  // dxdotdw[71]
    ddppS2_S4_S4dt_dflux_v27_v_26 = 0.075930144267274111;  // dxdotdw[72]
    ddS4dt_dflux_v28_v_27 = 0.15186028853454822;  // dxdotdw[73]
    ddpS2dt_dflux_v28_v_27 = 0.075930144267274111;  // dxdotdw[74]
    ddppS2_S4_S4dt_dflux_v28_v_27 = -0.075930144267274111;  // dxdotdw[75]
    ddS4dt_dflux_v29_v_28 = -0.15186028853454822;  // dxdotdw[76]
    ddppS3dt_dflux_v29_v_28 = -0.075930144267274111;  // dxdotdw[77]
    ddppS3_S4_S4dt_dflux_v29_v_28 = 0.075930144267274111;  // dxdotdw[78]
    ddS4dt_dflux_v30_v_29 = 0.15186028853454822;  // dxdotdw[79]
    ddpS3dt_dflux_v30_v_29 = 0.075930144267274111;  // dxdotdw[80]
    ddppS3_S4_S4dt_dflux_v30_v_29 = -0.075930144267274111;  // dxdotdw[81]
    ddS4dt_dflux_v31_v_30 = -0.075930144267274111;  // dxdotdw[82]
    ddppS2dt_dflux_v31_v_30 = -0.075930144267274111;  // dxdotdw[83]
    ddppS3dt_dflux_v31_v_30 = -0.075930144267274111;  // dxdotdw[84]
    ddppS2_ppS3_S4dt_dflux_v31_v_30 = 0.075930144267274111;  // dxdotdw[85]
    ddS4dt_dflux_v32_v_31 = 0.075930144267274111;  // dxdotdw[86]
    ddpS2dt_dflux_v32_v_31 = 0.075930144267274111;  // dxdotdw[87]
    ddppS3dt_dflux_v32_v_31 = 0.075930144267274111;  // dxdotdw[88]
    ddppS2_ppS3_S4dt_dflux_v32_v_31 = -0.075930144267274111;  // dxdotdw[89]
    ddS4dt_dflux_v33_v_32 = 0.075930144267274111;  // dxdotdw[90]
    ddpS3dt_dflux_v33_v_32 = 0.075930144267274111;  // dxdotdw[91]
    ddppS2dt_dflux_v33_v_32 = 0.075930144267274111;  // dxdotdw[92]
    ddppS2_ppS3_S4dt_dflux_v33_v_32 = -0.075930144267274111;  // dxdotdw[93]
    ddgeneAdt_dflux_v34_v_33 = 0.075930144267274111;  // dxdotdw[94]
    ddgeneAdt_dflux_v35_v_34 = -0.075930144267274111;  // dxdotdw[95]
    ddgeneBdt_dflux_v36_v_35 = 0.075930144267274111;  // dxdotdw[96]
    ddgeneBdt_dflux_v37_v_36 = -0.075930144267274111;  // dxdotdw[97]
    ddgeneCdt_dflux_v38_v_37 = 0.075930144267274111;  // dxdotdw[98]
    ddgeneCdt_dflux_v39_v_38 = -0.075930144267274111;  // dxdotdw[99]
    ddgeneDdt_dflux_v40_v_39 = 0.075930144267274111;  // dxdotdw[100]
    ddgeneDdt_dflux_v41_v_40 = -0.075930144267274111;  // dxdotdw[101]
    ddgeneEdt_dflux_v42_v_41 = 0.075930144267274111;  // dxdotdw[102]
    ddgeneEdt_dflux_v43_v_42 = -0.075930144267274111;  // dxdotdw[103]
    ddgeneFdt_dflux_v44_v_43 = 0.075930144267274111;  // dxdotdw[104]
    ddgeneFdt_dflux_v45_v_44 = -0.075930144267274111;  // dxdotdw[105]
    ddgeneGdt_dflux_v46_v_45 = 0.075930144267274111;  // dxdotdw[106]
    ddgeneGdt_dflux_v47_v_46 = -0.075930144267274111;  // dxdotdw[107]
    ddgeneHdt_dflux_v48_v_47 = 0.075930144267274111;  // dxdotdw[108]
    ddgeneHdt_dflux_v49_v_48 = -0.075930144267274111;  // dxdotdw[109]
    ddgeneIdt_dflux_v50_v_49 = 0.075930144267274111;  // dxdotdw[110]
    ddgeneIdt_dflux_v51_v_50 = -0.075930144267274111;  // dxdotdw[111]
    ddgeneJdt_dflux_v52_v_51 = 0.075930144267274111;  // dxdotdw[112]
    ddgeneJdt_dflux_v53_v_52 = -0.075930144267274111;  // dxdotdw[113]
    ddgeneKdt_dflux_v54_v_53 = 0.075930144267274111;  // dxdotdw[114]
    ddgeneKdt_dflux_v55_v_54 = -0.075930144267274111;  // dxdotdw[115]
    ddgeneLdt_dflux_v56_v_55 = 0.075930144267274111;  // dxdotdw[116]
    ddgeneLdt_dflux_v57_v_56 = -0.075930144267274111;  // dxdotdw[117]
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
