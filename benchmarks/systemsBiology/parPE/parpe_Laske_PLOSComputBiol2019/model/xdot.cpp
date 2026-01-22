#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void xdot_Laske_PLOSComputBiol2019(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dV_exdt = flux_V_att_Hi_Dissociation_virion + flux_V_att_Lo_Dissociation_virion - flux_V_ex_Virion_attachment_Hi - flux_V_ex_Virion_attachment_Lo;  // xdot[0]
    dV_enddt = flux_V_att_Hi_Endocytosis_virion + flux_V_att_Lo_Endocytosis_virion - flux_V_end_Degradation - flux_V_end_Fusion_and_release;  // xdot[1]
    dVp_cytdt = -flux_V_cyt_Import + 8*flux_V_end_Fusion_and_release;  // xdot[2]
    dVp_nucdt = flux_V_cyt_Import - flux_Vp_nuc_binding_of_P_M1 - flux_Vp_nuc_degradation + flux_Vp_nuc_synthesis;  // xdot[3]
    dR_Cdt = -flux_R_C_Degradation + flux_R_C_Synthesis - flux_R_C_binding_of_RdRp;  // xdot[4]
    dR_C_RdRpdt = -flux_C_P_synthesis - flux_R_C_RdRp_Degradation + flux_R_C_binding_of_RdRp;  // xdot[5]
    dR_Vdt = -flux_R_V_binding_of_RdRp - flux_R_V_degradation + flux_R_V_synthesis;  // xdot[6]
    dR_V_RdRpdt = -flux_R_V_RdRp_degradation + flux_R_V_binding_of_RdRp - flux_Vp_nuc_synthesis;  // xdot[7]
    dR_M1dt = -flux_R_M1_degradation + flux_R_M1_synthesis;  // xdot[8]
    dR_M2dt = -flux_R_M2_degradation + flux_R_M2_synthesis;  // xdot[9]
    dR_M3dt = -flux_R_M3_degradation + flux_R_M3_synthesis;  // xdot[10]
    dR_M4dt = -flux_R_M4_degradation + flux_R_M4_synthesis;  // xdot[11]
    dR_M5dt = -flux_R_M5_degradation + flux_R_M5_synthesis;  // xdot[12]
    dR_M6dt = -flux_R_M6_degradation + flux_R_M6_synthesis;  // xdot[13]
    dR_M7dt = -flux_R_M7_degradation + flux_R_M7_synthesis;  // xdot[14]
    dR_M8dt = -flux_R_M8_degradation + flux_R_M8_synthesis;  // xdot[15]
    dP_B1dt = flux_P_B1_synthesis - flux_P_RdRp_complex_formation;  // xdot[16]
    dP_B2dt = flux_P_B2_synthesis - flux_P_RdRp_complex_formation;  // xdot[17]
    dP_PAdt = flux_P_PA_synthesis - flux_P_RdRp_complex_formation;  // xdot[18]
    dP_HAdt = flux_P_HA_synthesis - 500*flux_Vp_cyt_M1_Virion_release;  // xdot[19]
    dP_NPdt = -71*flux_C_P_synthesis + flux_P_NP_synthesis - 433*flux_Vp_cyt_M1_Virion_release - 71*flux_Vp_nuc_synthesis;  // xdot[20]
    dP_NAdt = flux_P_NA_synthesis - 100*flux_Vp_cyt_M1_Virion_release;  // xdot[21]
    dP_M1dt = flux_P_M1_synthesis - 2932*flux_Vp_cyt_M1_Virion_release - 8.5*flux_Vp_nuc_binding_of_P_M1;  // xdot[22]
    dP_NEPdt = flux_P_NEP_synthesis - 157*flux_Vp_cyt_M1_Virion_release - flux_Vp_nuc_M1_export;  // xdot[23]
    dVp_nuc_M1dt = -flux_Vp_nuc_M1_degradation - flux_Vp_nuc_M1_export + flux_Vp_nuc_binding_of_P_M1;  // xdot[24]
    dVp_cyt_M1dt = -8*flux_Vp_cyt_M1_Virion_release - flux_Vp_cyt_M1_degradation + flux_Vp_nuc_M1_export;  // xdot[25]
    dV_reldt = flux_Vp_cyt_M1_Virion_release;  // xdot[26]
    dP_RdRpdt = flux_P_RdRp_complex_formation - flux_R_C_binding_of_RdRp - flux_R_V_binding_of_RdRp - 37*flux_Vp_cyt_M1_Virion_release;  // xdot[27]
    dCpdt = -flux_C_P_degradation + flux_C_P_synthesis;  // xdot[28]
    dV_att_Hidt = -flux_V_att_Hi_Dissociation_virion - flux_V_att_Hi_Endocytosis_virion + flux_V_ex_Virion_attachment_Hi;  // xdot[29]
    dV_att_Lodt = -flux_V_att_Lo_Dissociation_virion - flux_V_att_Lo_Endocytosis_virion + flux_V_ex_Virion_attachment_Lo;  // xdot[30]
    dP_M2dt = flux_P_M2_synthesis - 40*flux_Vp_cyt_M1_Virion_release;  // xdot[31]
    dB_att_Lodt = flux_V_att_Lo_Dissociation_virion + flux_V_att_Lo_Endocytosis_virion - flux_V_ex_Virion_attachment_Lo;  // xdot[32]
    dB_att_Hidt = flux_V_att_Hi_Dissociation_virion + flux_V_att_Hi_Endocytosis_virion - flux_V_ex_Virion_attachment_Hi;  // xdot[33]
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
