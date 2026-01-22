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

namespace amici {
namespace model_Raimundez_PCB2020 {

void w_Raimundez_PCB2020(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        AKT_total = AKTtotal__MKN1_fm*std::exp(HS746T*d_AKTtotal__MKN1_2_HS746T + d_AKTtotal__fm_2_hm*hunger_medium);  // w[0]
        FM = full_medium;  // w[3]
        HM = hunger_medium;  // w[4]
        MAPK_total = MAPKtotal__MKN1_fm*std::exp(HS746T*d_MAPKtotal__MKN1_2_HS746T + d_MAPKtotal__fm_2_hm*hunger_medium);  // w[5]
        MPI3K_total = MKN1*MPI3Ktotal__MKN1_fm*std::exp(d_MPI3Ktotal__fm_2_hm*hunger_medium);  // w[7]
        PI3K_total = HS746T*PI3Ktotal__HS746T_fm*std::exp(d_PI3Ktotal__fm_2_hm*hunger_medium);  // w[8]
        RAS_total = RAStotal__MKN1_fm*std::exp(HS746T*d_RAStotal__MKN1_2_HS746T + d_RAStotal__fm_2_hm*hunger_medium);  // w[9]
        ka_AKT = ka_AKT__MKN1*std::exp(HS746T*d_ka_AKT__MKN1_2_HS746T);  // w[10]
        ka_MAPK = ka_MAPK__MKN1*std::exp(HS746T*d_ka_MAPK__MKN1_2_HS746T);  // w[11]
        ka_RAS__basal = ka_RAS__basal__MKN1*std::exp(HS746T*d_ka_RAS__basal__MKN1_2_HS746T);  // w[12]
        ka_RAS__pEGFR_EGF_2 = ka_RAS__pEGFR_EGF_2__MKN1*std::exp(HS746T*d_ka_RAS__pEGFR_EGF_2__MKN1_2_HS746T);  // w[13]
        kdeg_membran = kdeg_membran__MKN1*std::exp(HS746T*d_kdeg_membran__MKN1_2_HS746T);  // w[14]
        kdeg_pEGFR_EGF_2_i = kdeg_pEGFR_EGF_2_i__MKN1*std::exp(HS746T*d_kdeg_pEGFR_EGF_2_i__MKN1_2_HS746T);  // w[15]
        kexp_pEGFR_EGF_2_i = kexp_pEGFR_EGF_2_i__MKN1*std::exp(HS746T*d_kexp_pEGFR_EGF_2_i__MKN1_2_HS746T);  // w[16]
        ki_AKT = ki_AKT__MKN1*std::exp(HS746T*d_ki_AKT__MKN1_2_HS746T);  // w[17]
        ki_PI3K = ki_PI3K__MKN1*std::exp(HS746T*d_ki_PI3K__MKN1_2_HS746T);  // w[18]
        ki_RAS = ki_RAS__MKN1*std::exp(HS746T*d_ki_RAS__MKN1_2_HS746T);  // w[19]
        kimp_pEGFR_EGF_2 = kimp_pEGFR_EGF_2__MKN1*std::exp(HS746T*d_kimp_pEGFR_EGF_2__MKN1_2_HS746T);  // w[20]
        ksyn_EGFR = ksyn_EGFR__MKN1_fm*relative_ksyn_EGFR*std::exp(HS746T*d_ksyn_EGFR__MKN1_2_HS746T + d_ksyn_EGFR__fm_2_hm*hunger_medium);  // w[21]
        ksyn_MMET = HS746T*ksyn_MMET__HS746T_fm*std::exp(d_ksyn_MMET__fm_2_hm*hunger_medium);  // w[22]
        flux_v1_v_0 = 1.0*ksyn_EGFR;  // w[23]
        flux_v15_v_14 = 1.0*ksyn_MMET;  // w[37]
    }

    // dynamic expressions
    CET = 6.8600000000000003*CET_level*Heaviside_0;  // w[1]
    EGF = 0.161*EGF_level*Heaviside_0;  // w[2]
    METinh = Heaviside_0*METinh_level;  // w[6]
    flux_v2_v_1 = 1.0*EGFR*kdeg_membran;  // w[24]
    flux_v3_v_2 = 1.0*CET*EGFR*kbin_EGFR_CET;  // w[25]
    flux_v4_v_3 = 1.0*EGFR_CET*KD_EGFR_CET*kbin_EGFR_CET;  // w[26]
    flux_v5_v_4 = 1.0*EGFR_CET*kdeg_membran;  // w[27]
    flux_v6_v_5 = 1.0*EGF*EGFR*kbin_EGFR_EGF;  // w[28]
    flux_v7_v_6 = 1.0*EGFR_EGF*KD_EGFR_EGF*kbin_EGFR_EGF;  // w[29]
    flux_v8_v_7 = 1.0*std::pow(EGFR_EGF, 2)*kdim_EGFR_EGF;  // w[30]
    flux_v9_v_8 = 1.0*EGFR_EGF_2*kpho_EGFR_EGF;  // w[31]
    flux_v10_v_9 = 1.0*kimp_pEGFR_EGF_2*pEGFR_EGF_2;  // w[32]
    flux_v11_v_10 = 1.0*kexp_pEGFR_EGF_2_i*pEGFR_EGF_2_i;  // w[33]
    flux_v12_v_11 = 1.0*EGFR_EGF*kdeg_membran;  // w[34]
    flux_v13_v_12 = 1.0*EGFR_EGF_2*kdeg_membran;  // w[35]
    flux_v14_v_13 = 1.0*kdeg_pEGFR_EGF_2_i*pEGFR_EGF_2_i;  // w[36]
    flux_v16_v_15 = 1.0*MMET*kdeg_membran;  // w[38]
    flux_v17_v_16 = 1.0*METinh*MMET*kdim_MMETinh;  // w[39]
    flux_v18_v_17 = 1.0*KD_METinh*MMET_METinh*kdim_MMETinh;  // w[40]
    flux_v19_v_18 = 1.0*MMET_METinh*kdeg_membran;  // w[41]
    flux_v20_v_19 = 1.0*std::pow(MMET, 2)*kdim_MMET*xi_kdim_MMET;  // w[42]
    flux_v21_v_20 = 1.0*MMET_2*kpho_MMET*xi_kpho_MMET;  // w[43]
    flux_v22_v_21 = 1.0*kimp_pEGFR_EGF_2*pMMET_2*xi_kimp_pMMET_2;  // w[44]
    flux_v23_v_22 = 1.0*kexp_pEGFR_EGF_2_i*pMMET_2_i*xi_kexp_pMMET_2_i;  // w[45]
    flux_v24_v_23 = 1.0*MMET_2*kdeg_membran;  // w[46]
    flux_v25_v_24 = 1.0*kdeg_pEGFR_EGF_2_i*pMMET_2_i*xi_kdeg_pMMET_2_i;  // w[47]
    flux_v26_v_25 = 1.0*MMET*MMET_METinh*kdim_MMET*xi_kdim_MMET;  // w[48]
    flux_v27_v_26 = 1.0*std::pow(MMET_METinh, 2)*kdim_MMET*xi_kdim_MMET;  // w[49]
    flux_v28_v_27 = 1.0*METinh*MMET_2*kdim_MMETinh;  // w[50]
    flux_v29_v_28 = 1.0*METinh*MMET_MMET_METinh*kdim_MMETinh;  // w[51]
    flux_v30_v_29 = 1.0*MMET_MMET_METinh*kdeg_membran;  // w[52]
    flux_v31_v_30 = 1.0*MMET_METinh_2*kdeg_membran;  // w[53]
    flux_v32_v_31 = 1.0*EGFR*MMET*kdim_MMET_EGFR*xi_kdim_MMET_EGFR;  // w[54]
    flux_v33_v_32 = 1.0*MMET_EGFR*kpho_MMET_EGFR*xi_kpho_MMET_EGFR;  // w[55]
    flux_v34_v_33 = 1.0*kimp_pEGFR_EGF_2*pMMET_pEGFR*xi_kimp_pMMET_pEGFR;  // w[56]
    flux_v35_v_34 = 1.0*kexp_pEGFR_EGF_2_i*pMMET_pEGFR_i*xi_kexp_pMMET_pEGFR_i;  // w[57]
    flux_v36_v_35 = 1.0*kdeg_pEGFR_EGF_2_i*pMMET_pEGFR_i*xi_kdeg_pMMET_pEGFR_i;  // w[58]
    flux_v37_v_36 = 1.0*MMET_EGFR*kdeg_membran;  // w[59]
    flux_v38_v_37 = 1.0*EGFR*MMET_METinh*kdim_MMET_EGFR*xi_kdim_MMET_EGFR;  // w[60]
    flux_v39_v_38 = -1.0*ka_RAS__basal*(RAS_GTP - RAS_total);  // w[61]
    flux_v40_v_39 = -1.0*ka_RAS__pEGFR_EGF_2*pEGFR_EGF_2*(RAS_GTP - RAS_total);  // w[62]
    flux_v41_v_40 = -1.0*ka_RAS__pEGFR_EGF_2*pEGFR_EGF_2_i*(RAS_GTP - RAS_total);  // w[63]
    flux_v42_v_41 = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_RAS_pMMET_pEGFR*(RAS_GTP - RAS_total);  // w[64]
    flux_v43_v_42 = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_RAS_pMMET_pEGFR*(RAS_GTP - RAS_total);  // w[65]
    flux_v44_v_43 = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_2*xi_ka_RAS_pMMET_2*(RAS_GTP - RAS_total);  // w[66]
    flux_v45_v_44 = -1.0*ka_RAS__pEGFR_EGF_2*pMMET_2_i*xi_ka_RAS_pMMET_2*(RAS_GTP - RAS_total);  // w[67]
    flux_v46_v_45 = 1.0*RAS_GTP*ki_RAS;  // w[68]
    flux_v47_v_46 = 1.0*RAS_GTP*ka_MAPK*(MAPK_total - pMAPK);  // w[69]
    flux_v48_v_47 = 1.0*ki_MAPK*pMAPK;  // w[70]
    flux_v49_v_48 = 1.0*ka_PI3K__basal*(PI3K_total - pPI3K);  // w[71]
    flux_v50_v_49 = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2*(PI3K_total - pPI3K);  // w[72]
    flux_v51_v_50 = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2_i*(PI3K_total - pPI3K);  // w[73]
    flux_v52_v_51 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR*(PI3K_total - pPI3K);  // w[74]
    flux_v53_v_52 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR*(PI3K_total - pPI3K);  // w[75]
    flux_v54_v_53 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*xi_ka_PI3K_pMMET_2*(PI3K_total - pPI3K);  // w[76]
    flux_v55_v_54 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*xi_ka_PI3K_pMMET_2*(PI3K_total - pPI3K);  // w[77]
    flux_v56_v_55 = 1.0*ki_PI3K*pPI3K;  // w[78]
    flux_v57_v_56 = 1.0*ka_PI3K__basal*(MPI3K_total - pMPI3K);  // w[79]
    flux_v58_v_57 = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2*(MPI3K_total - pMPI3K);  // w[80]
    flux_v59_v_58 = 1.0*ka_PI3K__pEGFR_EGF_2*pEGFR_EGF_2_i*(MPI3K_total - pMPI3K);  // w[81]
    flux_v60_v_59 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR*xi_ka_PI3K_pMMET_pEGFR*(MPI3K_total - pMPI3K);  // w[82]
    flux_v61_v_60 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_pEGFR_i*xi_ka_PI3K_pMMET_pEGFR*(MPI3K_total - pMPI3K);  // w[83]
    flux_v62_v_61 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2*xi_ka_PI3K_pMMET_2*(MPI3K_total - pMPI3K);  // w[84]
    flux_v63_v_62 = 1.0*ka_PI3K__pEGFR_EGF_2*pMMET_2_i*xi_ka_PI3K_pMMET_2*(MPI3K_total - pMPI3K);  // w[85]
    flux_v64_v_63 = 1.0*ki_PI3K*pMPI3K*xi_ki_MPI3K;  // w[86]
    flux_v65_v_64 = 1.0*ka_AKT*pPI3K*(AKT_total - pAKT);  // w[87]
    flux_v66_v_65 = 1.0*ka_AKT*pMPI3K*(AKT_total - pAKT);  // w[88]
    flux_v67_v_66 = 1.0*ki_AKT*pAKT;  // w[89]
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
