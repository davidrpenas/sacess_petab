#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 23> dwdx_colptrs_Raimundez_PCB2020_ = {
    0, 5, 7, 10, 12, 16, 21, 26, 29, 33, 38, 40, 44, 49, 54, 56, 57, 57, 66, 68, 77, 86, 89
};

void dwdx_colptrs_Raimundez_PCB2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Raimundez_PCB2020_));
}
} // namespace model_Raimundez_PCB2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 89> dwdx_rowvals_Raimundez_PCB2020_ = {
    24, 25, 28, 54, 60, 26, 27, 29, 30, 34, 31, 35, 32, 62, 72, 80, 33, 36, 63, 73, 81, 38, 39, 42, 48, 54, 43, 46, 50, 44, 66, 76, 84, 45, 47, 67, 77, 85, 55, 59, 56, 64, 74, 82, 57, 58, 65, 75, 83, 40, 41, 48, 49, 60, 51, 52, 53, 61, 62, 63, 64, 65, 66, 67, 68, 69, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 87, 79, 80, 81, 82, 83, 84, 85, 86, 88, 87, 88, 89
};

void dwdx_rowvals_Raimundez_PCB2020(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Raimundez_PCB2020_));
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
#include "dwdx.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void dwdx_Raimundez_PCB2020(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v2_v_1_dEGFR = 1.0*kdeg_membran;  // dwdx[0]
        dflux_v4_v_3_dEGFR_CET = 1.0*KD_EGFR_CET*kbin_EGFR_CET;  // dwdx[5]
        dflux_v5_v_4_dEGFR_CET = 1.0*kdeg_membran;  // dwdx[6]
        dflux_v7_v_6_dEGFR_EGF = 1.0*KD_EGFR_EGF*kbin_EGFR_EGF;  // dwdx[7]
        dflux_v12_v_11_dEGFR_EGF = 1.0*kdeg_membran;  // dwdx[9]
        dflux_v9_v_8_dEGFR_EGF_2 = 1.0*kpho_EGFR_EGF;  // dwdx[10]
        dflux_v13_v_12_dEGFR_EGF_2 = 1.0*kdeg_membran;  // dwdx[11]
        dflux_v10_v_9_dpEGFR_EGF_2 = 1.0*kimp_pEGFR_EGF_2;  // dwdx[12]
        dflux_v11_v_10_dpEGFR_EGF_2_i = 1.0*kexp_pEGFR_EGF_2_i;  // dwdx[16]
        dflux_v14_v_13_dpEGFR_EGF_2_i = 1.0*kdeg_pEGFR_EGF_2_i;  // dwdx[17]
        dflux_v16_v_15_dMMET = 1.0*kdeg_membran;  // dwdx[21]
        dflux_v21_v_20_dMMET_2 = 1.0*kpho_MMET*xi_kpho_MMET;  // dwdx[26]
        dflux_v24_v_23_dMMET_2 = 1.0*kdeg_membran;  // dwdx[27]
        dflux_v22_v_21_dpMMET_2 = 1.0*kimp_pEGFR_EGF_2*xi_kimp_pMMET_2;  // dwdx[29]
        dflux_v23_v_22_dpMMET_2_i = 1.0*kexp_pEGFR_EGF_2_i*xi_kexp_pMMET_2_i;  // dwdx[33]
        dflux_v25_v_24_dpMMET_2_i = 1.0*kdeg_pEGFR_EGF_2_i*xi_kdeg_pMMET_2_i;  // dwdx[34]
        dflux_v33_v_32_dMMET_EGFR = 1.0*kpho_MMET_EGFR*xi_kpho_MMET_EGFR;  // dwdx[38]
        dflux_v37_v_36_dMMET_EGFR = 1.0*kdeg_membran;  // dwdx[39]
        dflux_v34_v_33_dpMMET_pEGFR = 1.0*kimp_pEGFR_EGF_2*xi_kimp_pMMET_pEGFR;  // dwdx[40]
        dflux_v35_v_34_dpMMET_pEGFR_i = 1.0*kexp_pEGFR_EGF_2_i*xi_kexp_pMMET_pEGFR_i;  // dwdx[44]
        dflux_v36_v_35_dpMMET_pEGFR_i = 1.0*kdeg_pEGFR_EGF_2_i*xi_kdeg_pMMET_pEGFR_i;  // dwdx[45]
        dflux_v18_v_17_dMMET_METinh = 1.0*KD_METinh*kdim_MMETinh;  // dwdx[49]
        dflux_v19_v_18_dMMET_METinh = 1.0*kdeg_membran;  // dwdx[50]
        dflux_v30_v_29_dMMET_MMET_METinh = 1.0*kdeg_membran;  // dwdx[55]
        dflux_v31_v_30_dMMET_METinh_2 = 1.0*kdeg_membran;  // dwdx[56]
        dflux_v39_v_38_dRAS_GTP = -1.0*ka_RAS__basal;  // dwdx[57]
        dflux_v46_v_45_dRAS_GTP = 1.0*ki_RAS;  // dwdx[64]
        dflux_v48_v_47_dpMAPK = 1.0*ki_MAPK;  // dwdx[67]
        dflux_v49_v_48_dpPI3K = -1.0*ka_PI3K__basal;  // dwdx[68]
        dflux_v56_v_55_dpPI3K = 1.0*ki_PI3K;  // dwdx[75]
        dflux_v57_v_56_dpMPI3K = -1.0*ka_PI3K__basal;  // dwdx[77]
        dflux_v64_v_63_dpMPI3K = 1.0*ki_PI3K*xi_ki_MPI3K;  // dwdx[84]
        dflux_v67_v_66_dpAKT = 1.0*ki_AKT;  // dwdx[88]
    }

    // dynamic expressions
    dflux_v3_v_2_dEGFR = 1.0*CET*kbin_EGFR_CET;  // dwdx[1]
    dflux_v6_v_5_dEGFR = 1.0*EGF*kbin_EGFR_EGF;  // dwdx[2]
    dflux_v32_v_31_dEGFR = 1.0*MMET*kdim_MMET_EGFR*xi_kdim_MMET_EGFR;  // dwdx[3]
    dflux_v38_v_37_dEGFR = 1.0*MMET_METinh*kdim_MMET_EGFR*xi_kdim_MMET_EGFR;  // dwdx[4]
    dflux_v8_v_7_dEGFR_EGF = 2.0*EGFR_EGF*kdim_EGFR_EGF;  // dwdx[8]
    dflux_v40_v_39_dpEGFR_EGF_2 = -1.0*ka_RAS__pEGFR_EGF_2*(RAS_GTP - RAS_total);  // dwdx[13]
    dflux_v50_v_49_dpEGFR_EGF_2 = 1.0*ka_PI3K__pEGFR_EGF_2*(PI3K_total - pPI3K);  // dwdx[14]
    dflux_v58_v_57_dpEGFR_EGF_2 = 1.0*ka_PI3K__pEGFR_EGF_2*(MPI3K_total - pMPI3K);  // dwdx[15]
    dflux_v41_v_40_dpEGFR_EGF_2_i = -1.0*ka_RAS__pEGFR_EGF_2*(RAS_GTP - RAS_total);  // dwdx[18]
    dflux_v51_v_50_dpEGFR_EGF_2_i = 1.0*ka_PI3K__pEGFR_EGF_2*(PI3K_total - pPI3K);  // dwdx[19]
    dflux_v59_v_58_dpEGFR_EGF_2_i = 1.0*ka_PI3K__pEGFR_EGF_2*(MPI3K_total - pMPI3K);  // dwdx[20]
    dflux_v17_v_16_dMMET = 1.0*METinh*kdim_MMETinh;  // dwdx[22]
    dflux_v20_v_19_dMMET = 2.0*MMET*kdim_MMET*xi_kdim_MMET;  // dwdx[23]
    dflux_v26_v_25_dMMET = 1.0*MMET_METinh*kdim_MMET*xi_kdim_MMET;  // dwdx[24]
    dflux_v32_v_31_dMMET = 1.0*EGFR*kdim_MMET_EGFR*xi_kdim_MMET_EGFR;  // dwdx[25]
    dflux_v28_v_27_dMMET_2 = 1.0*METinh*kdim_MMETinh;  // dwdx[28]
    dflux_v44_v_43_dpMMET_2 = -1.0*ka_RAS__pEGFR_EGF_2*xi_ka_RAS_pMMET_2*(RAS_GTP - RAS_total);  // dwdx[30]
    dflux_v54_v_53_dpMMET_2 = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_2*(PI3K_total - pPI3K);  // dwdx[31]
    dflux_v62_v_61_dpMMET_2 = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_2*(MPI3K_total - pMPI3K);  // dwdx[32]
    dflux_v45_v_44_dpMMET_2_i = -1.0*ka_RAS__pEGFR_EGF_2*xi_ka_RAS_pMMET_2*(RAS_GTP - RAS_total);  // dwdx[35]
    dflux_v55_v_54_dpMMET_2_i = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_2*(PI3K_total - pPI3K);  // dwdx[36]
    dflux_v63_v_62_dpMMET_2_i = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_2*(MPI3K_total - pMPI3K);  // dwdx[37]
    dflux_v42_v_41_dpMMET_pEGFR = -1.0*ka_RAS__pEGFR_EGF_2*xi_ka_RAS_pMMET_pEGFR*(RAS_GTP - RAS_total);  // dwdx[41]
    dflux_v52_v_51_dpMMET_pEGFR = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_pEGFR*(PI3K_total - pPI3K);  // dwdx[42]
    dflux_v60_v_59_dpMMET_pEGFR = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_pEGFR*(MPI3K_total - pMPI3K);  // dwdx[43]
    dflux_v43_v_42_dpMMET_pEGFR_i = -1.0*ka_RAS__pEGFR_EGF_2*xi_ka_RAS_pMMET_pEGFR*(RAS_GTP - RAS_total);  // dwdx[46]
    dflux_v53_v_52_dpMMET_pEGFR_i = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_pEGFR*(PI3K_total - pPI3K);  // dwdx[47]
    dflux_v61_v_60_dpMMET_pEGFR_i = 1.0*ka_PI3K__pEGFR_EGF_2*xi_ka_PI3K_pMMET_pEGFR*(MPI3K_total - pMPI3K);  // dwdx[48]
    dflux_v26_v_25_dMMET_METinh = 1.0*MMET*kdim_MMET*xi_kdim_MMET;  // dwdx[51]
    dflux_v27_v_26_dMMET_METinh = 2.0*MMET_METinh*kdim_MMET*xi_kdim_MMET;  // dwdx[52]
    dflux_v38_v_37_dMMET_METinh = 1.0*EGFR*kdim_MMET_EGFR*xi_kdim_MMET_EGFR;  // dwdx[53]
    dflux_v29_v_28_dMMET_MMET_METinh = 1.0*METinh*kdim_MMETinh;  // dwdx[54]
    dflux_v40_v_39_dRAS_GTP = -1.0*ka_RAS__pEGFR_EGF_2*pEGFR_EGF_2;  // dwdx[58]
    dflux_v41_v_40_dRAS_GTP = -1.0*ka_RAS__pEGFR_EGF_2*pEGFR_EGF_2_i;  // dwdx[59]
    dflux_v42_v_41_dRAS_GTP = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_RAS_pMMET_pEGFR;  // dwdx[60]
    dflux_v43_v_42_dRAS_GTP = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_RAS_pMMET_pEGFR;  // dwdx[61]
    dflux_v44_v_43_dRAS_GTP = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_2*xi_ka_RAS_pMMET_2;  // dwdx[62]
    dflux_v45_v_44_dRAS_GTP = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_2_i*xi_ka_RAS_pMMET_2;  // dwdx[63]
    dflux_v47_v_46_dRAS_GTP = 1.0*ka_MAPK*(MAPK_total - pMAPK);  // dwdx[65]
    dflux_v47_v_46_dpMAPK = -1.0*RAS_GTP*ka_MAPK;  // dwdx[66]
    dflux_v50_v_49_dpPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2;  // dwdx[69]
    dflux_v51_v_50_dpPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2_i;  // dwdx[70]
    dflux_v52_v_51_dpPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR;  // dwdx[71]
    dflux_v53_v_52_dpPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR;  // dwdx[72]
    dflux_v54_v_53_dpPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*xi_ka_PI3K_pMMET_2;  // dwdx[73]
    dflux_v55_v_54_dpPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*xi_ka_PI3K_pMMET_2;  // dwdx[74]
    dflux_v65_v_64_dpPI3K = 1.0*ka_AKT*(AKT_total - pAKT);  // dwdx[76]
    dflux_v58_v_57_dpMPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2;  // dwdx[78]
    dflux_v59_v_58_dpMPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2_i;  // dwdx[79]
    dflux_v60_v_59_dpMPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR;  // dwdx[80]
    dflux_v61_v_60_dpMPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR;  // dwdx[81]
    dflux_v62_v_61_dpMPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*xi_ka_PI3K_pMMET_2;  // dwdx[82]
    dflux_v63_v_62_dpMPI3K = -1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*xi_ka_PI3K_pMMET_2;  // dwdx[83]
    dflux_v66_v_65_dpMPI3K = 1.0*ka_AKT*(AKT_total - pAKT);  // dwdx[85]
    dflux_v65_v_64_dpAKT = -1.0*ka_AKT*pPI3K;  // dwdx[86]
    dflux_v66_v_65_dpAKT = -1.0*ka_AKT*pMPI3K;  // dwdx[87]
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
