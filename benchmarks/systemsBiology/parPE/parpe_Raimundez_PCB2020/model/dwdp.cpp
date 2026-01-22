#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 248> dwdp_colptrs_Raimundez_PCB2020_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 34, 35, 36, 38, 40, 41, 42, 43, 46, 48, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 66, 70, 72, 74, 75, 76, 79, 81, 82, 83, 84, 85, 86, 87, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88
};

void dwdp_colptrs_Raimundez_PCB2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Raimundez_PCB2020_));
}
} // namespace model_Raimundez_PCB2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raimundez_PCB2020 {

static constexpr std::array<sunindextype, 88> dwdp_rowvals_Raimundez_PCB2020_ = {
    26, 29, 40, 0, 0, 5, 5, 7, 8, 9, 9, 14, 15, 16, 20, 21, 21, 22, 10, 11, 71, 79, 72, 73, 74, 75, 76, 77, 80, 81, 82, 83, 84, 85, 12, 13, 25, 26, 28, 29, 14, 15, 30, 42, 48, 49, 54, 60, 39, 40, 50, 51, 16, 17, 70, 18, 19, 20, 31, 43, 55, 22, 76, 77, 84, 85, 74, 75, 82, 83, 66, 67, 64, 65, 47, 58, 42, 48, 49, 54, 60, 45, 57, 86, 44, 56, 43, 55
};

void dwdp_rowvals_Raimundez_PCB2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Raimundez_PCB2020_));
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
#include "dwdp.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void dwdp_Raimundez_PCB2020(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dAKT_total_dd_AKTtotal__MKN1_2_HS746T = AKTtotal__MKN1_fm*HS746T*std::exp(HS746T*d_AKTtotal__MKN1_2_HS746T + d_AKTtotal__fm_2_hm*hunger_medium);  // dwdp[3]
        dAKT_total_dd_AKTtotal__fm_2_hm = AKTtotal__MKN1_fm*hunger_medium*std::exp(HS746T*d_AKTtotal__MKN1_2_HS746T + d_AKTtotal__fm_2_hm*hunger_medium);  // dwdp[4]
        dMAPK_total_dd_MAPKtotal__MKN1_2_HS746T = HS746T*MAPKtotal__MKN1_fm*std::exp(HS746T*d_MAPKtotal__MKN1_2_HS746T + d_MAPKtotal__fm_2_hm*hunger_medium);  // dwdp[5]
        dMAPK_total_dd_MAPKtotal__fm_2_hm = MAPKtotal__MKN1_fm*hunger_medium*std::exp(HS746T*d_MAPKtotal__MKN1_2_HS746T + d_MAPKtotal__fm_2_hm*hunger_medium);  // dwdp[6]
        dMPI3K_total_dd_MPI3Ktotal__fm_2_hm = MKN1*MPI3Ktotal__MKN1_fm*hunger_medium*std::exp(d_MPI3Ktotal__fm_2_hm*hunger_medium);  // dwdp[7]
        dPI3K_total_dd_PI3Ktotal__fm_2_hm = HS746T*PI3Ktotal__HS746T_fm*hunger_medium*std::exp(d_PI3Ktotal__fm_2_hm*hunger_medium);  // dwdp[8]
        dRAS_total_dd_RAStotal__MKN1_2_HS746T = HS746T*RAStotal__MKN1_fm*std::exp(HS746T*d_RAStotal__MKN1_2_HS746T + d_RAStotal__fm_2_hm*hunger_medium);  // dwdp[9]
        dRAS_total_dd_RAStotal__fm_2_hm = RAStotal__MKN1_fm*hunger_medium*std::exp(HS746T*d_RAStotal__MKN1_2_HS746T + d_RAStotal__fm_2_hm*hunger_medium);  // dwdp[10]
        dkdeg_membran_dd_kdeg_membran__MKN1_2_HS746T = HS746T*kdeg_membran__MKN1*std::exp(HS746T*d_kdeg_membran__MKN1_2_HS746T);  // dwdp[11]
        dkdeg_pEGFR_EGF_2_i_dd_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T = HS746T*kdeg_pEGFR_EGF_2_i__MKN1*std::exp(HS746T*d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T);  // dwdp[12]
        dkexp_pEGFR_EGF_2_i_dd_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T = HS746T*kexp_pEGFR_EGF_2_i__MKN1*std::exp(HS746T*d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T);  // dwdp[13]
        dkimp_pEGFR_EGF_2_dd_kimp_pEGFR_EGF_2__MKN1_2_HS746T = HS746T*kimp_pEGFR_EGF_2__MKN1*std::exp(HS746T*d_kimp_pEGFR_EGF_2__MKN1_2_HS746T);  // dwdp[14]
        dksyn_EGFR_dd_ksyn_EGFR__MKN1_2_HS746T = HS746T*ksyn_EGFR__MKN1_fm*relative_ksyn_EGFR*std::exp(HS746T*d_ksyn_EGFR__MKN1_2_HS746T + d_ksyn_EGFR__fm_2_hm*hunger_medium);  // dwdp[15]
        dksyn_EGFR_dd_ksyn_EGFR__fm_2_hm = hunger_medium*ksyn_EGFR__MKN1_fm*relative_ksyn_EGFR*std::exp(HS746T*d_ksyn_EGFR__MKN1_2_HS746T + d_ksyn_EGFR__fm_2_hm*hunger_medium);  // dwdp[16]
        dksyn_MMET_dd_ksyn_MMET__fm_2_hm = HS746T*hunger_medium*ksyn_MMET__HS746T_fm*std::exp(d_ksyn_MMET__fm_2_hm*hunger_medium);  // dwdp[17]
        dka_AKT_dka_AKT__MKN1 = std::exp(HS746T*d_ka_AKT__MKN1_2_HS746T);  // dwdp[18]
        dka_MAPK_dka_MAPK__MKN1 = std::exp(HS746T*d_ka_MAPK__MKN1_2_HS746T);  // dwdp[19]
        dka_RAS__basal_dka_RAS__basal__MKN1 = std::exp(HS746T*d_ka_RAS__basal__MKN1_2_HS746T);  // dwdp[34]
        dka_RAS__pEGFR_EGF_2_dka_RAS__pEGFR_EGF_2__MKN1 = std::exp(HS746T*d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T);  // dwdp[35]
        dkdeg_membran_dkdeg_membran__MKN1 = std::exp(HS746T*d_kdeg_membran__MKN1_2_HS746T);  // dwdp[40]
        dkdeg_pEGFR_EGF_2_i_dkdeg_pEGFR_EGF_2_i__MKN1 = std::exp(HS746T*d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T);  // dwdp[41]
        dkexp_pEGFR_EGF_2_i_dkexp_pEGFR_EGF_2_i__MKN1 = std::exp(HS746T*d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T);  // dwdp[52]
        dki_AKT_dki_AKT__MKN1 = std::exp(HS746T*d_ki_AKT__MKN1_2_HS746T);  // dwdp[53]
        dki_PI3K_dki_PI3K__MKN1 = std::exp(HS746T*d_ki_PI3K__MKN1_2_HS746T);  // dwdp[55]
        dki_RAS_dki_RAS__MKN1 = std::exp(HS746T*d_ki_RAS__MKN1_2_HS746T);  // dwdp[56]
        dkimp_pEGFR_EGF_2_dkimp_pEGFR_EGF_2__MKN1 = std::exp(HS746T*d_kimp_pEGFR_EGF_2__MKN1_2_HS746T);  // dwdp[57]
        dksyn_MMET_dksyn_MMET__HS746T_fm = HS746T*std::exp(d_ksyn_MMET__fm_2_hm*hunger_medium);  // dwdp[61]
    }

    // dynamic expressions
    dflux_v4_v_3_dKD_EGFR_CET = 1.0*EGFR_CET*kbin_EGFR_CET;  // dwdp[0]
    dflux_v7_v_6_dKD_EGFR_EGF = 1.0*EGFR_EGF*kbin_EGFR_EGF;  // dwdp[1]
    dflux_v18_v_17_dKD_METinh = 1.0*MMET_METinh*kdim_MMETinh;  // dwdp[2]
    dflux_v49_v_48_dka_PI3K__basal = 1.0*PI3K_total - 1.0*pPI3K;  // dwdp[20]
    dflux_v57_v_56_dka_PI3K__basal = 1.0*MPI3K_total - 1.0*pMPI3K;  // dwdp[21]
    dflux_v50_v_49_dka_PI3K__pEGFR_EGF_2 = 1.0*pEGFR_EGF_2*(PI3K_total - pPI3K);  // dwdp[22]
    dflux_v51_v_50_dka_PI3K__pEGFR_EGF_2 = 1.0*pEGFR_EGF_2_i*(PI3K_total - pPI3K);  // dwdp[23]
    dflux_v52_v_51_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR*(PI3K_total - pPI3K);  // dwdp[24]
    dflux_v53_v_52_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR*(PI3K_total - pPI3K);  // dwdp[25]
    dflux_v54_v_53_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_2*xi_ka_PI3K_pMMET_2*(PI3K_total - pPI3K);  // dwdp[26]
    dflux_v55_v_54_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_2_i*xi_ka_PI3K_pMMET_2*(PI3K_total - pPI3K);  // dwdp[27]
    dflux_v58_v_57_dka_PI3K__pEGFR_EGF_2 = 1.0*pEGFR_EGF_2*(MPI3K_total - pMPI3K);  // dwdp[28]
    dflux_v59_v_58_dka_PI3K__pEGFR_EGF_2 = 1.0*pEGFR_EGF_2_i*(MPI3K_total - pMPI3K);  // dwdp[29]
    dflux_v60_v_59_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR*(MPI3K_total - pMPI3K);  // dwdp[30]
    dflux_v61_v_60_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR*(MPI3K_total - pMPI3K);  // dwdp[31]
    dflux_v62_v_61_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_2*xi_ka_PI3K_pMMET_2*(MPI3K_total - pMPI3K);  // dwdp[32]
    dflux_v63_v_62_dka_PI3K__pEGFR_EGF_2 = 1.0*pMMET_2_i*xi_ka_PI3K_pMMET_2*(MPI3K_total - pMPI3K);  // dwdp[33]
    dflux_v3_v_2_dkbin_EGFR_CET = 1.0*CET*EGFR;  // dwdp[36]
    dflux_v4_v_3_dkbin_EGFR_CET = 1.0*EGFR_CET*KD_EGFR_CET;  // dwdp[37]
    dflux_v6_v_5_dkbin_EGFR_EGF = 1.0*EGF*EGFR;  // dwdp[38]
    dflux_v7_v_6_dkbin_EGFR_EGF = 1.0*EGFR_EGF*KD_EGFR_EGF;  // dwdp[39]
    dflux_v8_v_7_dkdim_EGFR_EGF = 1.0*std::pow(EGFR_EGF, 2);  // dwdp[42]
    dflux_v20_v_19_dkdim_MMET = 1.0*std::pow(MMET, 2)*xi_kdim_MMET;  // dwdp[43]
    dflux_v26_v_25_dkdim_MMET = 1.0*MMET*MMET_METinh*xi_kdim_MMET;  // dwdp[44]
    dflux_v27_v_26_dkdim_MMET = 1.0*std::pow(MMET_METinh, 2)*xi_kdim_MMET;  // dwdp[45]
    dflux_v32_v_31_dkdim_MMET_EGFR = 1.0*EGFR*MMET*xi_kdim_MMET_EGFR;  // dwdp[46]
    dflux_v38_v_37_dkdim_MMET_EGFR = 1.0*EGFR*MMET_METinh*xi_kdim_MMET_EGFR;  // dwdp[47]
    dflux_v17_v_16_dkdim_MMETinh = 1.0*METinh*MMET;  // dwdp[48]
    dflux_v18_v_17_dkdim_MMETinh = 1.0*KD_METinh*MMET_METinh;  // dwdp[49]
    dflux_v28_v_27_dkdim_MMETinh = 1.0*METinh*MMET_2;  // dwdp[50]
    dflux_v29_v_28_dkdim_MMETinh = 1.0*METinh*MMET_MMET_METinh;  // dwdp[51]
    dflux_v48_v_47_dki_MAPK = 1.0*pMAPK;  // dwdp[54]
    dflux_v9_v_8_dkpho_EGFR_EGF = 1.0*EGFR_EGF_2;  // dwdp[58]
    dflux_v21_v_20_dkpho_MMET = 1.0*MMET_2*xi_kpho_MMET;  // dwdp[59]
    dflux_v33_v_32_dkpho_MMET_EGFR = 1.0*MMET_EGFR*xi_kpho_MMET_EGFR;  // dwdp[60]
    dflux_v54_v_53_dxi_ka_PI3K_pMMET_2 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*(PI3K_total - pPI3K);  // dwdp[62]
    dflux_v55_v_54_dxi_ka_PI3K_pMMET_2 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*(PI3K_total - pPI3K);  // dwdp[63]
    dflux_v62_v_61_dxi_ka_PI3K_pMMET_2 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*(MPI3K_total - pMPI3K);  // dwdp[64]
    dflux_v63_v_62_dxi_ka_PI3K_pMMET_2 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*(MPI3K_total - pMPI3K);  // dwdp[65]
    dflux_v52_v_51_dxi_ka_PI3K_pMMET_pEGFR = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*(PI3K_total - pPI3K);  // dwdp[66]
    dflux_v53_v_52_dxi_ka_PI3K_pMMET_pEGFR = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*(PI3K_total - pPI3K);  // dwdp[67]
    dflux_v60_v_59_dxi_ka_PI3K_pMMET_pEGFR = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*(MPI3K_total - pMPI3K);  // dwdp[68]
    dflux_v61_v_60_dxi_ka_PI3K_pMMET_pEGFR = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*(MPI3K_total - pMPI3K);  // dwdp[69]
    dflux_v44_v_43_dxi_ka_RAS_pMMET_2 = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_2*(RAS_GTP - RAS_total);  // dwdp[70]
    dflux_v45_v_44_dxi_ka_RAS_pMMET_2 = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_2_i*(RAS_GTP - RAS_total);  // dwdp[71]
    dflux_v42_v_41_dxi_ka_RAS_pMMET_pEGFR = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR*(RAS_GTP - RAS_total);  // dwdp[72]
    dflux_v43_v_42_dxi_ka_RAS_pMMET_pEGFR = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR_i*(RAS_GTP - RAS_total);  // dwdp[73]
    dflux_v25_v_24_dxi_kdeg_pMMET_2_i = 1.0*kdeg_pEGFR_EGF_2_i*pMMET_2_i;  // dwdp[74]
    dflux_v36_v_35_dxi_kdeg_pMMET_pEGFR_i = 1.0*kdeg_pEGFR_EGF_2_i*pMMET_pEGFR_i;  // dwdp[75]
    dflux_v20_v_19_dxi_kdim_MMET = 1.0*std::pow(MMET, 2)*kdim_MMET;  // dwdp[76]
    dflux_v26_v_25_dxi_kdim_MMET = 1.0*MMET*MMET_METinh*kdim_MMET;  // dwdp[77]
    dflux_v27_v_26_dxi_kdim_MMET = 1.0*std::pow(MMET_METinh, 2)*kdim_MMET;  // dwdp[78]
    dflux_v32_v_31_dxi_kdim_MMET_EGFR = 1.0*EGFR*MMET*kdim_MMET_EGFR;  // dwdp[79]
    dflux_v38_v_37_dxi_kdim_MMET_EGFR = 1.0*EGFR*MMET_METinh*kdim_MMET_EGFR;  // dwdp[80]
    dflux_v23_v_22_dxi_kexp_pMMET_2_i = 1.0*kexp_pEGFR_EGF_2_i*pMMET_2_i;  // dwdp[81]
    dflux_v35_v_34_dxi_kexp_pMMET_pEGFR_i = 1.0*kexp_pEGFR_EGF_2_i*pMMET_pEGFR_i;  // dwdp[82]
    dflux_v64_v_63_dxi_ki_MPI3K = 1.0*ki_PI3K*pMPI3K;  // dwdp[83]
    dflux_v22_v_21_dxi_kimp_pMMET_2 = 1.0*kimp_pEGFR_EGF_2*pMMET_2;  // dwdp[84]
    dflux_v34_v_33_dxi_kimp_pMMET_pEGFR = 1.0*kimp_pEGFR_EGF_2*pMMET_pEGFR;  // dwdp[85]
    dflux_v21_v_20_dxi_kpho_MMET = 1.0*MMET_2*kpho_MMET;  // dwdp[86]
    dflux_v33_v_32_dxi_kpho_MMET_EGFR = 1.0*MMET_EGFR*kpho_MMET_EGFR;  // dwdp[87]
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
