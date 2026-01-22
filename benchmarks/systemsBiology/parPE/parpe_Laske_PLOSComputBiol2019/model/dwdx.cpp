#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 35> dwdx_colptrs_Laske_PLOSComputBiol2019_ = {
    0, 2, 6, 9, 22, 25, 28, 31, 33, 35, 37, 39, 41, 43, 45, 48, 50, 51, 52, 53, 54, 57, 58, 60, 62, 66, 70, 70, 73, 76, 79, 82, 83, 84, 85
};

void dwdx_colptrs_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Laske_PLOSComputBiol2019_));
}
} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 85> dwdx_rowvals_Laske_PLOSComputBiol2019_ = {
    56, 101, 27, 30, 59, 60, 27, 30, 61, 28, 30, 62, 72, 73, 74, 75, 76, 77, 78, 79, 80, 98, 29, 63, 64, 29, 65, 66, 30, 69, 70, 71, 106, 82, 89, 81, 90, 83, 91, 84, 92, 85, 93, 86, 94, 87, 95, 102, 88, 96, 97, 97, 97, 100, 66, 71, 100, 100, 98, 100, 99, 100, 28, 30, 99, 105, 27, 30, 100, 107, 64, 69, 100, 29, 67, 68, 30, 57, 58, 30, 103, 104, 100, 101, 56
};

void dwdx_rowvals_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Laske_PLOSComputBiol2019_));
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
#include "dwdx.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void dwdx_Laske_PLOSComputBiol2019(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dRNP_cyt_dV_end = 8;  // dwdx[2]
        dR_V_total_0_dV_end = 8;  // dwdx[3]
        dflux_V_end_Degradation_dV_end = 1.0*k_deg_end;  // dwdx[4]
        dflux_V_end_Fusion_and_release_dV_end = 1.0*k_fus;  // dwdx[5]
        dRNP_cyt_dVp_cyt = 1;  // dwdx[6]
        dR_V_total_0_dVp_cyt = 1;  // dwdx[7]
        dflux_V_cyt_Import_dVp_cyt = 1.0*k_imp;  // dwdx[8]
        dRNP_nuc_dVp_nuc = 1;  // dwdx[9]
        dR_V_total_0_dVp_nuc = 1;  // dwdx[10]
        dflux_R_C_Synthesis_dVp_nuc = 1.0*k_syn_R_C;  // dwdx[11]
        dflux_Vp_nuc_degradation_dVp_nuc = 1.0*k_deg_Rnp;  // dwdx[12]
        dflux_R_M1_synthesis_dVp_nuc = 1.0*k_syn_R_M1_rel;  // dwdx[13]
        dflux_R_M2_synthesis_dVp_nuc = 1.0*k_syn_R_M2_rel;  // dwdx[14]
        dflux_R_M3_synthesis_dVp_nuc = 1.0*k_syn_R_M3_rel;  // dwdx[15]
        dflux_R_M4_synthesis_dVp_nuc = 1.0*k_syn_R_M4_rel;  // dwdx[16]
        dflux_R_M5_synthesis_dVp_nuc = 1.0*k_syn_R_M5_rel;  // dwdx[17]
        dflux_R_M6_synthesis_dVp_nuc = 1.0*k_syn_R_M6_rel;  // dwdx[18]
        dflux_R_M7_synthesis_dVp_nuc = 1.0*k_syn_R_M7_rel;  // dwdx[19]
        dflux_R_M8_synthesis_dVp_nuc = 1.0*k_syn_R_M8_rel;  // dwdx[20]
        dR_C_total_0_dR_C = 1;  // dwdx[22]
        dflux_R_C_Degradation_dR_C = 1.0*k_deg_R;  // dwdx[23]
        dR_C_total_0_dR_C_RdRp = 1;  // dwdx[25]
        dflux_R_C_RdRp_Degradation_dR_C_RdRp = 1.0*k_deg_R_RdRp;  // dwdx[26]
        dR_V_total_0_dR_V = 1;  // dwdx[28]
        dflux_R_V_degradation_dR_V = 1.0*k_deg_R;  // dwdx[30]
        dflux_R_V_RdRp_degradation_dR_V_RdRp = 1.0*k_deg_R_RdRp;  // dwdx[32]
        dflux_P_B2_synthesis_dR_M1 = 1.0*k_syn_P_rel;  // dwdx[33]
        dflux_R_M1_degradation_dR_M1 = 1.0*k_deg_R_M;  // dwdx[34]
        dflux_P_B1_synthesis_dR_M2 = 1.0*k_syn_P_rel;  // dwdx[35]
        dflux_R_M2_degradation_dR_M2 = 1.0*k_deg_R_M;  // dwdx[36]
        dflux_P_PA_synthesis_dR_M3 = 1.0*k_syn_P_rel;  // dwdx[37]
        dflux_R_M3_degradation_dR_M3 = 1.0*k_deg_R_M;  // dwdx[38]
        dflux_P_HA_synthesis_dR_M4 = 1.0*k_syn_P_rel;  // dwdx[39]
        dflux_R_M4_degradation_dR_M4 = 1.0*k_deg_R_M;  // dwdx[40]
        dflux_P_NP_synthesis_dR_M5 = 1.0*k_syn_P_rel;  // dwdx[41]
        dflux_R_M5_degradation_dR_M5 = 1.0*k_deg_R_M;  // dwdx[42]
        dflux_P_NA_synthesis_dR_M6 = 1.0*k_syn_P_rel;  // dwdx[43]
        dflux_R_M6_degradation_dR_M6 = 1.0*k_deg_R_M;  // dwdx[44]
        dflux_P_M1_synthesis_dR_M7 = 1.0*k_syn_P_M1_rel;  // dwdx[45]
        dflux_R_M7_degradation_dR_M7 = 1.0*k_deg_R_M;  // dwdx[46]
        dflux_P_M2_synthesis_dR_M7 = 1.0*k_syn_P_M2_rel;  // dwdx[47]
        dflux_P_NEP_synthesis_dR_M8 = 1.0*k_syn_P_NEP_rel;  // dwdx[48]
        dflux_R_M8_degradation_dR_M8 = 1.0*k_deg_R_M;  // dwdx[49]
        dRNP_nuc_dVp_nuc_M1 = 1;  // dwdx[62]
        dR_V_total_0_dVp_nuc_M1 = 1;  // dwdx[63]
        dflux_Vp_nuc_M1_degradation_dVp_nuc_M1 = 1.0*k_deg_Rnp;  // dwdx[65]
        dRNP_cyt_dVp_cyt_M1 = 1;  // dwdx[66]
        dR_V_total_0_dVp_cyt_M1 = 1;  // dwdx[67]
        dflux_Vp_cyt_M1_degradation_dVp_cyt_M1 = 1.0*k_deg_Rnp;  // dwdx[69]
        dR_C_total_0_dCp = 1;  // dwdx[73]
        dflux_R_V_synthesis_dCp = 1.0*k_syn_R_V;  // dwdx[74]
        dflux_C_P_degradation_dCp = 1.0*k_deg_Rnp;  // dwdx[75]
        dR_V_total_0_dV_att_Hi = 8;  // dwdx[76]
        dflux_V_att_Hi_Dissociation_virion_dV_att_Hi = 1.0*k_dis_Hi;  // dwdx[77]
        dflux_V_att_Hi_Endocytosis_virion_dV_att_Hi = 1.0*k_end;  // dwdx[78]
        dR_V_total_0_dV_att_Lo = 8;  // dwdx[79]
        dflux_V_att_Lo_Dissociation_virion_dV_att_Lo = 1.0*k_dis_Lo;  // dwdx[80]
        dflux_V_att_Lo_Endocytosis_virion_dV_att_Lo = 1.0*k_end;  // dwdx[81]
    }

    // dynamic expressions
    dflux_V_ex_Virion_attachment_Hi_dV_ex = 1.0*B_att_Hi*k_att_Hi;  // dwdx[0]
    dflux_V_ex_Virion_attachment_Lo_dV_ex = 1.0*B_att_Lo*k_att_Lo;  // dwdx[1]
    dflux_Vp_nuc_binding_of_P_M1_dVp_nuc = 1.0*P_M1*k_bind_M1;  // dwdx[21]
    dflux_R_C_binding_of_RdRp_dR_C = 1.0*P_RdRp*k_bind_RdRp;  // dwdx[24]
    dflux_C_P_synthesis_dR_C_RdRp = 1.0*P_NP*k_bind_NP;  // dwdx[27]
    dflux_R_V_binding_of_RdRp_dR_V = 1.0*P_RdRp*k_bind_RdRp;  // dwdx[29]
    dflux_Vp_nuc_synthesis_dR_V_RdRp = 1.0*P_NP*k_bind_NP;  // dwdx[31]
    dflux_P_RdRp_complex_formation_dP_B1 = 1.0*P_B2*P_PA*k_RdRp;  // dwdx[50]
    dflux_P_RdRp_complex_formation_dP_B2 = 1.0*P_B1*P_PA*k_RdRp;  // dwdx[51]
    dflux_P_RdRp_complex_formation_dP_PA = 1.0*P_B1*P_B2*k_RdRp;  // dwdx[52]
    dflux_Vp_cyt_M1_Virion_release_dP_HA = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*std::pow(KmC + P_HA, 2)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP)) + 1.0*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[53]
    dflux_C_P_synthesis_dP_NP = 1.0*R_C_RdRp*k_bind_NP;  // dwdx[54]
    dflux_Vp_nuc_synthesis_dP_NP = 1.0*R_V_RdRp*k_bind_NP;  // dwdx[55]
    dflux_Vp_cyt_M1_Virion_release_dP_NP = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*std::pow(KmD + P_NP, 2)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP)) + 1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[56]
    dflux_Vp_cyt_M1_Virion_release_dP_NA = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*std::pow(KmE + P_NA, 2)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP)) + 1.0*P_HA*P_M1*P_M2*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[57]
    dflux_Vp_nuc_binding_of_P_M1_dP_M1 = 1.0*Vp_nuc*k_bind_M1;  // dwdx[58]
    dflux_Vp_cyt_M1_Virion_release_dP_M1 = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*std::pow(KmF + P_M1, 2)*(KmG + P_M2)*(KmH + P_NEP)) + 1.0*P_HA*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[59]
    dflux_Vp_nuc_M1_export_dP_NEP = 1.0*Vp_nuc_M1*k_exp_Vp_nuc_M1;  // dwdx[60]
    dflux_Vp_cyt_M1_Virion_release_dP_NEP = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*std::pow(KmH + P_NEP, 2)) + 1.0*P_HA*P_M1*P_M2*P_NA*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[61]
    dflux_Vp_nuc_M1_export_dVp_nuc_M1 = 1.0*P_NEP*k_exp_Vp_nuc_M1;  // dwdx[64]
    dflux_Vp_cyt_M1_Virion_release_dVp_cyt_M1 = 1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[68]
    dflux_R_C_binding_of_RdRp_dP_RdRp = 1.0*R_C*k_bind_RdRp;  // dwdx[70]
    dflux_R_V_binding_of_RdRp_dP_RdRp = 1.0*R_V*k_bind_RdRp;  // dwdx[71]
    dflux_Vp_cyt_M1_Virion_release_dP_RdRp = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/(std::pow(KmB + P_RdRp, 2)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP)) + 1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[72]
    dflux_Vp_cyt_M1_Virion_release_dP_M2 = -1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*std::pow(KmG + P_M2, 2)*(KmH + P_NEP)) + 1.0*P_HA*P_M1*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1*k_rel/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdx[82]
    dflux_V_ex_Virion_attachment_Lo_dB_att_Lo = 1.0*V_ex*k_att_Lo;  // dwdx[83]
    dflux_V_ex_Virion_attachment_Hi_dB_att_Hi = 1.0*V_ex*k_att_Hi;  // dwdx[84]
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
