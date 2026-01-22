#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 109> dxdotdw_colptrs_Laske_PLOSComputBiol2019_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 9, 10, 12, 14, 15, 16, 19, 20, 23, 24, 25, 28, 29, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 61, 64, 67, 76, 79, 80, 83, 86, 87, 88, 89
};

void dxdotdw_colptrs_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Laske_PLOSComputBiol2019_));
}
} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 89> dxdotdw_rowvals_Laske_PLOSComputBiol2019_ = {
    0, 29, 33, 0, 29, 33, 1, 29, 33, 1, 1, 2, 2, 3, 4, 4, 4, 5, 27, 5, 5, 20, 28, 6, 28, 6, 7, 27, 6, 3, 7, 20, 3, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 27, 3, 22, 24, 23, 24, 25, 19, 20, 21, 22, 23, 25, 26, 27, 31, 0, 30, 32, 31, 0, 30, 32, 1, 30, 32, 24, 7, 25
};

void dxdotdw_rowvals_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Laske_PLOSComputBiol2019_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void dxdotdw_Laske_PLOSComputBiol2019(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddV_exdt_dflux_V_ex_Virion_attachment_Hi = -1;  // dxdotdw[0]
    ddV_att_Hidt_dflux_V_ex_Virion_attachment_Hi = 1;  // dxdotdw[1]
    ddB_att_Hidt_dflux_V_ex_Virion_attachment_Hi = -1;  // dxdotdw[2]
    ddV_exdt_dflux_V_att_Hi_Dissociation_virion = 1;  // dxdotdw[3]
    ddV_att_Hidt_dflux_V_att_Hi_Dissociation_virion = -1;  // dxdotdw[4]
    ddB_att_Hidt_dflux_V_att_Hi_Dissociation_virion = 1;  // dxdotdw[5]
    ddV_enddt_dflux_V_att_Hi_Endocytosis_virion = 1;  // dxdotdw[6]
    ddV_att_Hidt_dflux_V_att_Hi_Endocytosis_virion = -1;  // dxdotdw[7]
    ddB_att_Hidt_dflux_V_att_Hi_Endocytosis_virion = 1;  // dxdotdw[8]
    ddV_enddt_dflux_V_end_Degradation = -1;  // dxdotdw[9]
    ddV_enddt_dflux_V_end_Fusion_and_release = -1;  // dxdotdw[10]
    ddVp_cytdt_dflux_V_end_Fusion_and_release = 8;  // dxdotdw[11]
    ddVp_cytdt_dflux_V_cyt_Import = -1;  // dxdotdw[12]
    ddVp_nucdt_dflux_V_cyt_Import = 1;  // dxdotdw[13]
    ddR_Cdt_dflux_R_C_Synthesis = 1;  // dxdotdw[14]
    ddR_Cdt_dflux_R_C_Degradation = -1;  // dxdotdw[15]
    ddR_Cdt_dflux_R_C_binding_of_RdRp = -1;  // dxdotdw[16]
    ddR_C_RdRpdt_dflux_R_C_binding_of_RdRp = 1;  // dxdotdw[17]
    ddP_RdRpdt_dflux_R_C_binding_of_RdRp = -1;  // dxdotdw[18]
    ddR_C_RdRpdt_dflux_R_C_RdRp_Degradation = -1;  // dxdotdw[19]
    ddR_C_RdRpdt_dflux_C_P_synthesis = -1;  // dxdotdw[20]
    ddP_NPdt_dflux_C_P_synthesis = -71;  // dxdotdw[21]
    ddCpdt_dflux_C_P_synthesis = 1;  // dxdotdw[22]
    ddR_Vdt_dflux_R_V_synthesis = 1;  // dxdotdw[23]
    ddCpdt_dflux_C_P_degradation = -1;  // dxdotdw[24]
    ddR_Vdt_dflux_R_V_binding_of_RdRp = -1;  // dxdotdw[25]
    ddR_V_RdRpdt_dflux_R_V_binding_of_RdRp = 1;  // dxdotdw[26]
    ddP_RdRpdt_dflux_R_V_binding_of_RdRp = -1;  // dxdotdw[27]
    ddR_Vdt_dflux_R_V_degradation = -1;  // dxdotdw[28]
    ddVp_nucdt_dflux_Vp_nuc_synthesis = 1;  // dxdotdw[29]
    ddR_V_RdRpdt_dflux_Vp_nuc_synthesis = -1;  // dxdotdw[30]
    ddP_NPdt_dflux_Vp_nuc_synthesis = -71;  // dxdotdw[31]
    ddVp_nucdt_dflux_Vp_nuc_degradation = -1;  // dxdotdw[32]
    ddR_M1dt_dflux_R_M1_synthesis = 1;  // dxdotdw[33]
    ddR_M2dt_dflux_R_M2_synthesis = 1;  // dxdotdw[34]
    ddR_M3dt_dflux_R_M3_synthesis = 1;  // dxdotdw[35]
    ddR_M4dt_dflux_R_M4_synthesis = 1;  // dxdotdw[36]
    ddR_M5dt_dflux_R_M5_synthesis = 1;  // dxdotdw[37]
    ddR_M6dt_dflux_R_M6_synthesis = 1;  // dxdotdw[38]
    ddR_M7dt_dflux_R_M7_synthesis = 1;  // dxdotdw[39]
    ddR_M8dt_dflux_R_M8_synthesis = 1;  // dxdotdw[40]
    ddP_B1dt_dflux_P_B1_synthesis = 1;  // dxdotdw[41]
    ddP_B2dt_dflux_P_B2_synthesis = 1;  // dxdotdw[42]
    ddP_PAdt_dflux_P_PA_synthesis = 1;  // dxdotdw[43]
    ddP_HAdt_dflux_P_HA_synthesis = 1;  // dxdotdw[44]
    ddP_NPdt_dflux_P_NP_synthesis = 1;  // dxdotdw[45]
    ddP_NAdt_dflux_P_NA_synthesis = 1;  // dxdotdw[46]
    ddP_M1dt_dflux_P_M1_synthesis = 1;  // dxdotdw[47]
    ddP_NEPdt_dflux_P_NEP_synthesis = 1;  // dxdotdw[48]
    ddR_M1dt_dflux_R_M1_degradation = -1;  // dxdotdw[49]
    ddR_M2dt_dflux_R_M2_degradation = -1;  // dxdotdw[50]
    ddR_M3dt_dflux_R_M3_degradation = -1;  // dxdotdw[51]
    ddR_M4dt_dflux_R_M4_degradation = -1;  // dxdotdw[52]
    ddR_M5dt_dflux_R_M5_degradation = -1;  // dxdotdw[53]
    ddR_M6dt_dflux_R_M6_degradation = -1;  // dxdotdw[54]
    ddR_M7dt_dflux_R_M7_degradation = -1;  // dxdotdw[55]
    ddR_M8dt_dflux_R_M8_degradation = -1;  // dxdotdw[56]
    ddP_B1dt_dflux_P_RdRp_complex_formation = -1;  // dxdotdw[57]
    ddP_B2dt_dflux_P_RdRp_complex_formation = -1;  // dxdotdw[58]
    ddP_PAdt_dflux_P_RdRp_complex_formation = -1;  // dxdotdw[59]
    ddP_RdRpdt_dflux_P_RdRp_complex_formation = 1;  // dxdotdw[60]
    ddVp_nucdt_dflux_Vp_nuc_binding_of_P_M1 = -1;  // dxdotdw[61]
    ddP_M1dt_dflux_Vp_nuc_binding_of_P_M1 = -8.5;  // dxdotdw[62]
    ddVp_nuc_M1dt_dflux_Vp_nuc_binding_of_P_M1 = 1;  // dxdotdw[63]
    ddP_NEPdt_dflux_Vp_nuc_M1_export = -1;  // dxdotdw[64]
    ddVp_nuc_M1dt_dflux_Vp_nuc_M1_export = -1;  // dxdotdw[65]
    ddVp_cyt_M1dt_dflux_Vp_nuc_M1_export = 1;  // dxdotdw[66]
    ddP_HAdt_dflux_Vp_cyt_M1_Virion_release = -500;  // dxdotdw[67]
    ddP_NPdt_dflux_Vp_cyt_M1_Virion_release = -433;  // dxdotdw[68]
    ddP_NAdt_dflux_Vp_cyt_M1_Virion_release = -100;  // dxdotdw[69]
    ddP_M1dt_dflux_Vp_cyt_M1_Virion_release = -2932;  // dxdotdw[70]
    ddP_NEPdt_dflux_Vp_cyt_M1_Virion_release = -157;  // dxdotdw[71]
    ddVp_cyt_M1dt_dflux_Vp_cyt_M1_Virion_release = -8;  // dxdotdw[72]
    ddV_reldt_dflux_Vp_cyt_M1_Virion_release = 1;  // dxdotdw[73]
    ddP_RdRpdt_dflux_Vp_cyt_M1_Virion_release = -37;  // dxdotdw[74]
    ddP_M2dt_dflux_Vp_cyt_M1_Virion_release = -40;  // dxdotdw[75]
    ddV_exdt_dflux_V_ex_Virion_attachment_Lo = -1;  // dxdotdw[76]
    ddV_att_Lodt_dflux_V_ex_Virion_attachment_Lo = 1;  // dxdotdw[77]
    ddB_att_Lodt_dflux_V_ex_Virion_attachment_Lo = -1;  // dxdotdw[78]
    ddP_M2dt_dflux_P_M2_synthesis = 1;  // dxdotdw[79]
    ddV_exdt_dflux_V_att_Lo_Dissociation_virion = 1;  // dxdotdw[80]
    ddV_att_Lodt_dflux_V_att_Lo_Dissociation_virion = -1;  // dxdotdw[81]
    ddB_att_Lodt_dflux_V_att_Lo_Dissociation_virion = 1;  // dxdotdw[82]
    ddV_enddt_dflux_V_att_Lo_Endocytosis_virion = 1;  // dxdotdw[83]
    ddV_att_Lodt_dflux_V_att_Lo_Endocytosis_virion = -1;  // dxdotdw[84]
    ddB_att_Lodt_dflux_V_att_Lo_Endocytosis_virion = 1;  // dxdotdw[85]
    ddVp_nuc_M1dt_dflux_Vp_nuc_M1_degradation = -1;  // dxdotdw[86]
    ddR_V_RdRpdt_dflux_R_V_RdRp_degradation = -1;  // dxdotdw[87]
    ddVp_cyt_M1dt_dflux_Vp_cyt_M1_degradation = -1;  // dxdotdw[88]
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
