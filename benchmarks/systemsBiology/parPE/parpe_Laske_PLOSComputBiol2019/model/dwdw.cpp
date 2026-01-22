#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 109> dwdw_colptrs_Laske_PLOSComputBiol2019_ = {
    0, 1, 2, 3, 4, 5, 6, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 31, 35, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 51, 52, 53, 54, 55, 56, 57, 58, 58, 58, 59, 60, 61, 62, 63, 64, 70, 71, 72, 73, 74, 75, 76, 77, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78
};

void dwdw_colptrs_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Laske_PLOSComputBiol2019_));
}
} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 78> dwdw_rowvals_Laske_PLOSComputBiol2019_ = {
    32, 33, 36, 35, 38, 34, 32, 33, 34, 35, 36, 37, 38, 44, 45, 46, 37, 42, 43, 41, 42, 43, 41, 48, 49, 50, 51, 52, 53, 54, 55, 44, 45, 46, 47, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 31, 31, 39, 40, 100, 100, 100, 100, 100, 100, 100, 59, 57, 103, 87, 102, 88, 81, 82, 83, 84, 85, 86, 73, 74, 75, 76, 77, 78, 79, 80
};

void dwdw_rowvals_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Laske_PLOSComputBiol2019_));
}
} // namespace model_Laske_PLOSComputBiol2019
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
#include "dwdw.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void dwdw_Laske_PLOSComputBiol2019(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dKmB_dModelValue_101 = ModelValue_107;  // dwdw[0]
        dKmC_dModelValue_102 = ModelValue_107;  // dwdw[1]
        dKmF_dModelValue_103 = ModelValue_107;  // dwdw[2]
        dKmE_dModelValue_104 = ModelValue_107;  // dwdw[3]
        dKmH_dModelValue_105 = ModelValue_107;  // dwdw[4]
        dKmD_dModelValue_106 = ModelValue_107;  // dwdw[5]
        dKmB_dModelValue_107 = ModelValue_101;  // dwdw[6]
        dKmC_dModelValue_107 = ModelValue_102;  // dwdw[7]
        dKmD_dModelValue_107 = ModelValue_106;  // dwdw[8]
        dKmE_dModelValue_107 = ModelValue_104;  // dwdw[9]
        dKmF_dModelValue_107 = ModelValue_103;  // dwdw[10]
        dKmG_dModelValue_107 = ModelValue_113;  // dwdw[11]
        dKmH_dModelValue_107 = ModelValue_105;  // dwdw[12]
        dk_syn_P_M1_rel_dModelValue_108 = -ModelValue_80/ModelValue_82;  // dwdw[13]
        dk_syn_P_M2_rel_dModelValue_108 = ModelValue_80/ModelValue_82;  // dwdw[14]
        dk_syn_P_NEP_rel_dModelValue_111 = ModelValue_80/ModelValue_82;  // dwdw[15]
        dKmG_dModelValue_113 = ModelValue_107;  // dwdw[16]
        dk_dis_Hi_dModelValue_114 = -ModelValue_63/std::pow(ModelValue_114, 2);  // dwdw[17]
        dk_dis_Lo_dModelValue_115 = -ModelValue_64/std::pow(ModelValue_115, 2);  // dwdw[18]
        dk_deg_end_dModelValue_116 = -ModelValue_69/ModelValue_116 - ModelValue_69*(1 - ModelValue_116)/std::pow(ModelValue_116, 2);  // dwdw[19]
        dk_dis_Hi_dModelValue_63 = 1.0/ModelValue_114;  // dwdw[20]
        dk_dis_Lo_dModelValue_64 = 1.0/ModelValue_115;  // dwdw[21]
        dk_deg_end_dModelValue_69 = (1 - ModelValue_116)/ModelValue_116;  // dwdw[22]
        dk_syn_R_M1_rel_dModelValue_79 = (1.0/8.0)/ModelValue_84;  // dwdw[23]
        dk_syn_R_M2_rel_dModelValue_79 = (1.0/8.0)/ModelValue_85;  // dwdw[24]
        dk_syn_R_M3_rel_dModelValue_79 = (1.0/8.0)/ModelValue_86;  // dwdw[25]
        dk_syn_R_M4_rel_dModelValue_79 = (1.0/8.0)/ModelValue_87;  // dwdw[26]
        dk_syn_R_M5_rel_dModelValue_79 = (1.0/8.0)/ModelValue_88;  // dwdw[27]
        dk_syn_R_M6_rel_dModelValue_79 = (1.0/8.0)/ModelValue_89;  // dwdw[28]
        dk_syn_R_M7_rel_dModelValue_79 = (1.0/8.0)/ModelValue_90;  // dwdw[29]
        dk_syn_R_M8_rel_dModelValue_79 = (1.0/8.0)/ModelValue_91;  // dwdw[30]
        dk_syn_P_M1_rel_dModelValue_80 = (1 - ModelValue_108)/ModelValue_82;  // dwdw[31]
        dk_syn_P_M2_rel_dModelValue_80 = ModelValue_108/ModelValue_82;  // dwdw[32]
        dk_syn_P_NEP_rel_dModelValue_80 = ModelValue_111/ModelValue_82;  // dwdw[33]
        dk_syn_P_rel_dModelValue_80 = 1.0/ModelValue_82;  // dwdw[34]
        dk_syn_P_M1_rel_dModelValue_82 = -ModelValue_80*(1 - ModelValue_108)/std::pow(ModelValue_82, 2);  // dwdw[35]
        dk_syn_P_M2_rel_dModelValue_82 = -ModelValue_108*ModelValue_80/std::pow(ModelValue_82, 2);  // dwdw[36]
        dk_syn_P_NEP_rel_dModelValue_82 = -ModelValue_111*ModelValue_80/std::pow(ModelValue_82, 2);  // dwdw[37]
        dk_syn_P_rel_dModelValue_82 = -ModelValue_80/std::pow(ModelValue_82, 2);  // dwdw[38]
        dk_syn_R_M1_rel_dModelValue_84 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_84, 2);  // dwdw[39]
        dk_syn_R_M2_rel_dModelValue_85 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_85, 2);  // dwdw[40]
        dk_syn_R_M3_rel_dModelValue_86 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_86, 2);  // dwdw[41]
        dk_syn_R_M4_rel_dModelValue_87 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_87, 2);  // dwdw[42]
        dk_syn_R_M5_rel_dModelValue_88 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_88, 2);  // dwdw[43]
        dk_syn_R_M6_rel_dModelValue_89 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_89, 2);  // dwdw[44]
        dk_syn_R_M7_rel_dModelValue_90 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_90, 2);  // dwdw[45]
        dk_syn_R_M8_rel_dModelValue_91 = -1.0/8.0*ModelValue_79/std::pow(ModelValue_91, 2);  // dwdw[46]
        dR_C_seg_tot_dR_C_total_0 = 1.0/8.0;  // dwdw[49]
        dR_V_seg_tot_dR_V_total_0 = 1.0/8.0;  // dwdw[50]
    }

    // dynamic expressions
    dF_rnp_nuc_dRNP_cyt = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);  // dwdw[47]
    dF_rnp_nuc_dRNP_nuc = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);  // dwdw[48]
    dflux_Vp_cyt_M1_Virion_release_dKmB = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/(std::pow(KmB + P_RdRp, 2)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdw[51]
    dflux_Vp_cyt_M1_Virion_release_dKmC = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*std::pow(KmC + P_HA, 2)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdw[52]
    dflux_Vp_cyt_M1_Virion_release_dKmD = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*std::pow(KmD + P_NP, 2)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdw[53]
    dflux_Vp_cyt_M1_Virion_release_dKmE = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*std::pow(KmE + P_NA, 2)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdw[54]
    dflux_Vp_cyt_M1_Virion_release_dKmF = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*std::pow(KmF + P_M1, 2)*(KmG + P_M2)*(KmH + P_NEP));  // dwdw[55]
    dflux_Vp_cyt_M1_Virion_release_dKmG = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*std::pow(KmG + P_M2, 2)*(KmH + P_NEP));  // dwdw[56]
    dflux_Vp_cyt_M1_Virion_release_dKmH = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*std::pow(KmH + P_NEP, 2));  // dwdw[57]
    dflux_V_end_Degradation_dk_deg_end = 1.0*V_end;  // dwdw[58]
    dflux_V_att_Hi_Dissociation_virion_dk_dis_Hi = 1.0*V_att_Hi;  // dwdw[59]
    dflux_V_att_Lo_Dissociation_virion_dk_dis_Lo = 1.0*V_att_Lo;  // dwdw[60]
    dflux_P_M1_synthesis_dk_syn_P_M1_rel = 1.0*R_M7;  // dwdw[61]
    dflux_P_M2_synthesis_dk_syn_P_M2_rel = 1.0*R_M7;  // dwdw[62]
    dflux_P_NEP_synthesis_dk_syn_P_NEP_rel = 1.0*R_M8;  // dwdw[63]
    dflux_P_B1_synthesis_dk_syn_P_rel = 1.0*R_M2;  // dwdw[64]
    dflux_P_B2_synthesis_dk_syn_P_rel = 1.0*R_M1;  // dwdw[65]
    dflux_P_PA_synthesis_dk_syn_P_rel = 1.0*R_M3;  // dwdw[66]
    dflux_P_HA_synthesis_dk_syn_P_rel = 1.0*R_M4;  // dwdw[67]
    dflux_P_NP_synthesis_dk_syn_P_rel = 1.0*R_M5;  // dwdw[68]
    dflux_P_NA_synthesis_dk_syn_P_rel = 1.0*R_M6;  // dwdw[69]
    dflux_R_M1_synthesis_dk_syn_R_M1_rel = 1.0*Vp_nuc;  // dwdw[70]
    dflux_R_M2_synthesis_dk_syn_R_M2_rel = 1.0*Vp_nuc;  // dwdw[71]
    dflux_R_M3_synthesis_dk_syn_R_M3_rel = 1.0*Vp_nuc;  // dwdw[72]
    dflux_R_M4_synthesis_dk_syn_R_M4_rel = 1.0*Vp_nuc;  // dwdw[73]
    dflux_R_M5_synthesis_dk_syn_R_M5_rel = 1.0*Vp_nuc;  // dwdw[74]
    dflux_R_M6_synthesis_dk_syn_R_M6_rel = 1.0*Vp_nuc;  // dwdw[75]
    dflux_R_M7_synthesis_dk_syn_R_M7_rel = 1.0*Vp_nuc;  // dwdw[76]
    dflux_R_M8_synthesis_dk_syn_R_M8_rel = 1.0*Vp_nuc;  // dwdw[77]
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
