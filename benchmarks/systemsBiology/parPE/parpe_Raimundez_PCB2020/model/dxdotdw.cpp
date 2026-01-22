#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 91> dxdotdw_colptrs_Raimundez_PCB2020_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4, 6, 7, 9, 11, 13, 15, 17, 19, 20, 21, 22, 23, 24, 26, 28, 29, 31, 33, 35, 37, 38, 39, 42, 44, 46, 48, 49, 50, 53, 55, 57, 60, 61, 62, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94
};

void dxdotdw_colptrs_Raimundez_PCB2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Raimundez_PCB2020_));
}
} // namespace model_Raimundez_PCB2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 94> dxdotdw_rowvals_Raimundez_PCB2020_ = {
    0, 0, 0, 1, 0, 1, 1, 0, 2, 0, 2, 2, 3, 3, 4, 4, 5, 0, 5, 2, 3, 5, 6, 6, 6, 13, 6, 13, 13, 6, 7, 7, 8, 8, 9, 6, 9, 7, 9, 6, 13, 14, 13, 15, 7, 14, 14, 15, 14, 15, 0, 6, 10, 10, 11, 11, 12, 0, 6, 12, 12, 10, 0, 13, 16, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21
};

void dxdotdw_rowvals_Raimundez_PCB2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Raimundez_PCB2020_));
}
} // namespace model_Raimundez_PCB2020
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dxdotdw.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void dxdotdw_Raimundez_PCB2020(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddEGFRdt_dflux_v1_v_0 = 1;  // dxdotdw[0]
    ddEGFRdt_dflux_v2_v_1 = -1;  // dxdotdw[1]
    ddEGFRdt_dflux_v3_v_2 = -1;  // dxdotdw[2]
    ddEGFR_CETdt_dflux_v3_v_2 = 1;  // dxdotdw[3]
    ddEGFRdt_dflux_v4_v_3 = 1;  // dxdotdw[4]
    ddEGFR_CETdt_dflux_v4_v_3 = -1;  // dxdotdw[5]
    ddEGFR_CETdt_dflux_v5_v_4 = -1;  // dxdotdw[6]
    ddEGFRdt_dflux_v6_v_5 = -1;  // dxdotdw[7]
    ddEGFR_EGFdt_dflux_v6_v_5 = 1;  // dxdotdw[8]
    ddEGFRdt_dflux_v7_v_6 = 1;  // dxdotdw[9]
    ddEGFR_EGFdt_dflux_v7_v_6 = -1;  // dxdotdw[10]
    ddEGFR_EGFdt_dflux_v8_v_7 = -2;  // dxdotdw[11]
    ddEGFR_EGF_2dt_dflux_v8_v_7 = 1;  // dxdotdw[12]
    ddEGFR_EGF_2dt_dflux_v9_v_8 = -1;  // dxdotdw[13]
    ddpEGFR_EGF_2dt_dflux_v9_v_8 = 1;  // dxdotdw[14]
    ddpEGFR_EGF_2dt_dflux_v10_v_9 = -1;  // dxdotdw[15]
    ddpEGFR_EGF_2_idt_dflux_v10_v_9 = 1;  // dxdotdw[16]
    ddEGFRdt_dflux_v11_v_10 = 2;  // dxdotdw[17]
    ddpEGFR_EGF_2_idt_dflux_v11_v_10 = -1;  // dxdotdw[18]
    ddEGFR_EGFdt_dflux_v12_v_11 = -1;  // dxdotdw[19]
    ddEGFR_EGF_2dt_dflux_v13_v_12 = -1;  // dxdotdw[20]
    ddpEGFR_EGF_2_idt_dflux_v14_v_13 = -1;  // dxdotdw[21]
    ddMMETdt_dflux_v15_v_14 = 1;  // dxdotdw[22]
    ddMMETdt_dflux_v16_v_15 = -1;  // dxdotdw[23]
    ddMMETdt_dflux_v17_v_16 = -1;  // dxdotdw[24]
    ddMMET_METinhdt_dflux_v17_v_16 = 1;  // dxdotdw[25]
    ddMMETdt_dflux_v18_v_17 = 1;  // dxdotdw[26]
    ddMMET_METinhdt_dflux_v18_v_17 = -1;  // dxdotdw[27]
    ddMMET_METinhdt_dflux_v19_v_18 = -1;  // dxdotdw[28]
    ddMMETdt_dflux_v20_v_19 = -2;  // dxdotdw[29]
    ddMMET_2dt_dflux_v20_v_19 = 1;  // dxdotdw[30]
    ddMMET_2dt_dflux_v21_v_20 = -1;  // dxdotdw[31]
    ddpMMET_2dt_dflux_v21_v_20 = 1;  // dxdotdw[32]
    ddpMMET_2dt_dflux_v22_v_21 = -1;  // dxdotdw[33]
    ddpMMET_2_idt_dflux_v22_v_21 = 1;  // dxdotdw[34]
    ddMMETdt_dflux_v23_v_22 = 2;  // dxdotdw[35]
    ddpMMET_2_idt_dflux_v23_v_22 = -1;  // dxdotdw[36]
    ddMMET_2dt_dflux_v24_v_23 = -1;  // dxdotdw[37]
    ddpMMET_2_idt_dflux_v25_v_24 = -1;  // dxdotdw[38]
    ddMMETdt_dflux_v26_v_25 = -1;  // dxdotdw[39]
    ddMMET_METinhdt_dflux_v26_v_25 = -1;  // dxdotdw[40]
    ddMMET_MMET_METinhdt_dflux_v26_v_25 = 1;  // dxdotdw[41]
    ddMMET_METinhdt_dflux_v27_v_26 = -2;  // dxdotdw[42]
    ddMMET_METinh_2dt_dflux_v27_v_26 = 1;  // dxdotdw[43]
    ddMMET_2dt_dflux_v28_v_27 = -1;  // dxdotdw[44]
    ddMMET_MMET_METinhdt_dflux_v28_v_27 = 1;  // dxdotdw[45]
    ddMMET_MMET_METinhdt_dflux_v29_v_28 = -1;  // dxdotdw[46]
    ddMMET_METinh_2dt_dflux_v29_v_28 = 1;  // dxdotdw[47]
    ddMMET_MMET_METinhdt_dflux_v30_v_29 = -1;  // dxdotdw[48]
    ddMMET_METinh_2dt_dflux_v31_v_30 = -1;  // dxdotdw[49]
    ddEGFRdt_dflux_v32_v_31 = -1;  // dxdotdw[50]
    ddMMETdt_dflux_v32_v_31 = -1;  // dxdotdw[51]
    ddMMET_EGFRdt_dflux_v32_v_31 = 1;  // dxdotdw[52]
    ddMMET_EGFRdt_dflux_v33_v_32 = -1;  // dxdotdw[53]
    ddpMMET_pEGFRdt_dflux_v33_v_32 = 1;  // dxdotdw[54]
    ddpMMET_pEGFRdt_dflux_v34_v_33 = -1;  // dxdotdw[55]
    ddpMMET_pEGFR_idt_dflux_v34_v_33 = 1;  // dxdotdw[56]
    ddEGFRdt_dflux_v35_v_34 = 1;  // dxdotdw[57]
    ddMMETdt_dflux_v35_v_34 = 1;  // dxdotdw[58]
    ddpMMET_pEGFR_idt_dflux_v35_v_34 = -1;  // dxdotdw[59]
    ddpMMET_pEGFR_idt_dflux_v36_v_35 = -1;  // dxdotdw[60]
    ddMMET_EGFRdt_dflux_v37_v_36 = -1;  // dxdotdw[61]
    ddEGFRdt_dflux_v38_v_37 = -1;  // dxdotdw[62]
    ddMMET_METinhdt_dflux_v38_v_37 = -1;  // dxdotdw[63]
    ddEGFR_MMET_METinhdt_dflux_v38_v_37 = 1;  // dxdotdw[64]
    ddRAS_GTPdt_dflux_v39_v_38 = 1;  // dxdotdw[65]
    ddRAS_GTPdt_dflux_v40_v_39 = 1;  // dxdotdw[66]
    ddRAS_GTPdt_dflux_v41_v_40 = 1;  // dxdotdw[67]
    ddRAS_GTPdt_dflux_v42_v_41 = 1;  // dxdotdw[68]
    ddRAS_GTPdt_dflux_v43_v_42 = 1;  // dxdotdw[69]
    ddRAS_GTPdt_dflux_v44_v_43 = 1;  // dxdotdw[70]
    ddRAS_GTPdt_dflux_v45_v_44 = 1;  // dxdotdw[71]
    ddRAS_GTPdt_dflux_v46_v_45 = -1;  // dxdotdw[72]
    ddpMAPKdt_dflux_v47_v_46 = 1;  // dxdotdw[73]
    ddpMAPKdt_dflux_v48_v_47 = -1;  // dxdotdw[74]
    ddpPI3Kdt_dflux_v49_v_48 = 1;  // dxdotdw[75]
    ddpPI3Kdt_dflux_v50_v_49 = 1;  // dxdotdw[76]
    ddpPI3Kdt_dflux_v51_v_50 = 1;  // dxdotdw[77]
    ddpPI3Kdt_dflux_v52_v_51 = 1;  // dxdotdw[78]
    ddpPI3Kdt_dflux_v53_v_52 = 1;  // dxdotdw[79]
    ddpPI3Kdt_dflux_v54_v_53 = 1;  // dxdotdw[80]
    ddpPI3Kdt_dflux_v55_v_54 = 1;  // dxdotdw[81]
    ddpPI3Kdt_dflux_v56_v_55 = -1;  // dxdotdw[82]
    ddpMPI3Kdt_dflux_v57_v_56 = 1;  // dxdotdw[83]
    ddpMPI3Kdt_dflux_v58_v_57 = 1;  // dxdotdw[84]
    ddpMPI3Kdt_dflux_v59_v_58 = 1;  // dxdotdw[85]
    ddpMPI3Kdt_dflux_v60_v_59 = 1;  // dxdotdw[86]
    ddpMPI3Kdt_dflux_v61_v_60 = 1;  // dxdotdw[87]
    ddpMPI3Kdt_dflux_v62_v_61 = 1;  // dxdotdw[88]
    ddpMPI3Kdt_dflux_v63_v_62 = 1;  // dxdotdw[89]
    ddpMPI3Kdt_dflux_v64_v_63 = -1;  // dxdotdw[90]
    ddpAKTdt_dflux_v65_v_64 = 1;  // dxdotdw[91]
    ddpAKTdt_dflux_v66_v_65 = 1;  // dxdotdw[92]
    ddpAKTdt_dflux_v67_v_66 = -1;  // dxdotdw[93]
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
