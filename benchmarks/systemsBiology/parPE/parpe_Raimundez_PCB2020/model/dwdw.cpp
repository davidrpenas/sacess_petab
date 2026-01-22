#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 91> dwdw_colptrs_Raimundez_PCB2020_ = {
    0, 2, 3, 4, 4, 4, 5, 8, 15, 22, 29, 31, 32, 33, 39, 49, 52, 55, 56, 58, 59, 62, 63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
};

void dwdw_colptrs_Raimundez_PCB2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Raimundez_PCB2020_));
}
} // namespace model_Raimundez_PCB2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 64> dwdw_rowvals_Raimundez_PCB2020_ = {
    87, 88, 25, 28, 69, 39, 50, 51, 79, 80, 81, 82, 83, 84, 85, 71, 72, 73, 74, 75, 76, 77, 61, 62, 63, 64, 65, 66, 67, 87, 88, 69, 61, 62, 63, 64, 65, 66, 67, 24, 27, 34, 35, 38, 41, 46, 52, 53, 59, 36, 47, 58, 33, 45, 57, 89, 78, 86, 68, 32, 44, 56, 23, 37
};

void dwdw_rowvals_Raimundez_PCB2020(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Raimundez_PCB2020_));
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
#include "dwdw.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void dwdw_Raimundez_PCB2020(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v57_v_56_dMPI3K_total = 1.0*ka_PI3K__basal;  // dwdw[8]
        dflux_v49_v_48_dPI3K_total = 1.0*ka_PI3K__basal;  // dwdw[15]
        dflux_v39_v_38_dRAS_total = 1.0*ka_RAS__basal;  // dwdw[22]
        dflux_v1_v_0_dksyn_EGFR = 1.0;  // dwdw[62]
        dflux_v15_v_14_dksyn_MMET = 1.0;  // dwdw[63]
    }

    // dynamic expressions
    dflux_v65_v_64_dAKT_total = 1.0*ka_AKT*pPI3K;  // dwdw[0]
    dflux_v66_v_65_dAKT_total = 1.0*ka_AKT*pMPI3K;  // dwdw[1]
    dflux_v3_v_2_dCET = 1.0*EGFR*kbin_EGFR_CET;  // dwdw[2]
    dflux_v6_v_5_dEGF = 1.0*EGFR*kbin_EGFR_EGF;  // dwdw[3]
    dflux_v47_v_46_dMAPK_total = 1.0*RAS_GTP*ka_MAPK;  // dwdw[4]
    dflux_v17_v_16_dMETinh = 1.0*MMET*kdim_MMETinh;  // dwdw[5]
    dflux_v28_v_27_dMETinh = 1.0*MMET_2*kdim_MMETinh;  // dwdw[6]
    dflux_v29_v_28_dMETinh = 1.0*MMET_MMET_METinh*kdim_MMETinh;  // dwdw[7]
    dflux_v58_v_57_dMPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2;  // dwdw[9]
    dflux_v59_v_58_dMPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2_i;  // dwdw[10]
    dflux_v60_v_59_dMPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR;  // dwdw[11]
    dflux_v61_v_60_dMPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR;  // dwdw[12]
    dflux_v62_v_61_dMPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*xi_ka_PI3K_pMMET_2;  // dwdw[13]
    dflux_v63_v_62_dMPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*xi_ka_PI3K_pMMET_2;  // dwdw[14]
    dflux_v50_v_49_dPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2;  // dwdw[16]
    dflux_v51_v_50_dPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2_i;  // dwdw[17]
    dflux_v52_v_51_dPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR;  // dwdw[18]
    dflux_v53_v_52_dPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR;  // dwdw[19]
    dflux_v54_v_53_dPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*xi_ka_PI3K_pMMET_2;  // dwdw[20]
    dflux_v55_v_54_dPI3K_total = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*xi_ka_PI3K_pMMET_2;  // dwdw[21]
    dflux_v40_v_39_dRAS_total = 1.0*ka_RAS__pEGFR_EGF_2*pEGFR_EGF_2;  // dwdw[23]
    dflux_v41_v_40_dRAS_total = 1.0*ka_RAS__pEGFR_EGF_2*pEGFR_EGF_2_i;  // dwdw[24]
    dflux_v42_v_41_dRAS_total = 1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_RAS_pMMET_pEGFR;  // dwdw[25]
    dflux_v43_v_42_dRAS_total = 1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_RAS_pMMET_pEGFR;  // dwdw[26]
    dflux_v44_v_43_dRAS_total = 1.0*ka_RAS__pEGFR_EGF_2*pMMET_2*xi_ka_RAS_pMMET_2;  // dwdw[27]
    dflux_v45_v_44_dRAS_total = 1.0*ka_RAS__pEGFR_EGF_2*pMMET_2_i*xi_ka_RAS_pMMET_2;  // dwdw[28]
    dflux_v65_v_64_dka_AKT = 1.0*pPI3K*(AKT_total - pAKT);  // dwdw[29]
    dflux_v66_v_65_dka_AKT = 1.0*pMPI3K*(AKT_total - pAKT);  // dwdw[30]
    dflux_v47_v_46_dka_MAPK = 1.0*RAS_GTP*(MAPK_total - pMAPK);  // dwdw[31]
    dflux_v39_v_38_dka_RAS__basal = -1.0*RAS_GTP + 1.0*RAS_total;  // dwdw[32]
    dflux_v40_v_39_dka_RAS__pEGFR_EGF_2 = -1.0*pEGFR_EGF_2*(RAS_GTP - RAS_total);  // dwdw[33]
    dflux_v41_v_40_dka_RAS__pEGFR_EGF_2 = -1.0*pEGFR_EGF_2_i*(RAS_GTP - RAS_total);  // dwdw[34]
    dflux_v42_v_41_dka_RAS__pEGFR_EGF_2 = -1.0*pMMET_pEGFR*xi_ka_RAS_pMMET_pEGFR*(RAS_GTP - RAS_total);  // dwdw[35]
    dflux_v43_v_42_dka_RAS__pEGFR_EGF_2 = -1.0*pMMET_pEGFR_i*xi_ka_RAS_pMMET_pEGFR*(RAS_GTP - RAS_total);  // dwdw[36]
    dflux_v44_v_43_dka_RAS__pEGFR_EGF_2 = -1.0*pMMET_2*xi_ka_RAS_pMMET_2*(RAS_GTP - RAS_total);  // dwdw[37]
    dflux_v45_v_44_dka_RAS__pEGFR_EGF_2 = -1.0*pMMET_2_i*xi_ka_RAS_pMMET_2*(RAS_GTP - RAS_total);  // dwdw[38]
    dflux_v2_v_1_dkdeg_membran = 1.0*EGFR;  // dwdw[39]
    dflux_v5_v_4_dkdeg_membran = 1.0*EGFR_CET;  // dwdw[40]
    dflux_v12_v_11_dkdeg_membran = 1.0*EGFR_EGF;  // dwdw[41]
    dflux_v13_v_12_dkdeg_membran = 1.0*EGFR_EGF_2;  // dwdw[42]
    dflux_v16_v_15_dkdeg_membran = 1.0*MMET;  // dwdw[43]
    dflux_v19_v_18_dkdeg_membran = 1.0*MMET_METinh;  // dwdw[44]
    dflux_v24_v_23_dkdeg_membran = 1.0*MMET_2;  // dwdw[45]
    dflux_v30_v_29_dkdeg_membran = 1.0*MMET_MMET_METinh;  // dwdw[46]
    dflux_v31_v_30_dkdeg_membran = 1.0*MMET_METinh_2;  // dwdw[47]
    dflux_v37_v_36_dkdeg_membran = 1.0*MMET_EGFR;  // dwdw[48]
    dflux_v14_v_13_dkdeg_pEGFR_EGF_2_i = 1.0*pEGFR_EGF_2_i;  // dwdw[49]
    dflux_v25_v_24_dkdeg_pEGFR_EGF_2_i = 1.0*pMMET_2_i*xi_kdeg_pMMET_2_i;  // dwdw[50]
    dflux_v36_v_35_dkdeg_pEGFR_EGF_2_i = 1.0*pMMET_pEGFR_i*xi_kdeg_pMMET_pEGFR_i;  // dwdw[51]
    dflux_v11_v_10_dkexp_pEGFR_EGF_2_i = 1.0*pEGFR_EGF_2_i;  // dwdw[52]
    dflux_v23_v_22_dkexp_pEGFR_EGF_2_i = 1.0*pMMET_2_i*xi_kexp_pMMET_2_i;  // dwdw[53]
    dflux_v35_v_34_dkexp_pEGFR_EGF_2_i = 1.0*pMMET_pEGFR_i*xi_kexp_pMMET_pEGFR_i;  // dwdw[54]
    dflux_v67_v_66_dki_AKT = 1.0*pAKT;  // dwdw[55]
    dflux_v56_v_55_dki_PI3K = 1.0*pPI3K;  // dwdw[56]
    dflux_v64_v_63_dki_PI3K = 1.0*pMPI3K*xi_ki_MPI3K;  // dwdw[57]
    dflux_v46_v_45_dki_RAS = 1.0*RAS_GTP;  // dwdw[58]
    dflux_v10_v_9_dkimp_pEGFR_EGF_2 = 1.0*pEGFR_EGF_2;  // dwdw[59]
    dflux_v22_v_21_dkimp_pEGFR_EGF_2 = 1.0*pMMET_2*xi_kimp_pMMET_2;  // dwdw[60]
    dflux_v34_v_33_dkimp_pEGFR_EGF_2 = 1.0*pMMET_pEGFR*xi_kimp_pMMET_pEGFR;  // dwdw[61]
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
