#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<sunindextype, 40> dwdp_colptrs_Raia_CancerResearch2011_ = {
    0, 1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 15, 16, 17, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19
};

void dwdp_colptrs_Raia_CancerResearch2011(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Raia_CancerResearch2011_));
}
} // namespace model_Raia_CancerResearch2011
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<sunindextype, 19> dwdp_rowvals_Raia_CancerResearch2011_ = {
    15, 11, 4, 5, 4, 5, 0, 1, 3, 2, 13, 14, 13, 12, 9, 8, 7, 6, 10
};

void dwdp_rowvals_Raia_CancerResearch2011(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Raia_CancerResearch2011_));
}
} // namespace model_Raia_CancerResearch2011
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
#include "dwdp.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void dwdp_Raia_CancerResearch2011(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v16_v_15_dCD274mRNA_production = 10.0*pSTAT5;  // dwdp[0]
    dflux_v12_v_11_dDecoyR_binding = 226.5*DecoyR*il13_level;  // dwdp[1]
    dflux_v5_v_4_dJAK2_p_inhibition = -100.0*IL13_Rec*JAK2*JAK2_phosphorylation*SOCS3/std::pow(JAK2_p_inhibition*SOCS3 + 1, 2);  // dwdp[2]
    dflux_v6_v_5_dJAK2_p_inhibition = -100.0*JAK2*JAK2_phosphorylation*SOCS3*p_IL13_Rec/std::pow(JAK2_p_inhibition*SOCS3 + 1, 2);  // dwdp[3]
    dflux_v5_v_4_dJAK2_phosphorylation = 100.0*IL13_Rec*JAK2/(JAK2_p_inhibition*SOCS3 + 1);  // dwdp[4]
    dflux_v6_v_5_dJAK2_phosphorylation = 100.0*JAK2*p_IL13_Rec/(JAK2_p_inhibition*SOCS3 + 1);  // dwdp[5]
    dflux_v1_v_0_dKon_IL13Rec = 226.5*Rec*il13_level;  // dwdp[6]
    dflux_v2_v_1_dRec_intern = 100.0*Rec;  // dwdp[7]
    dflux_v4_v_3_dRec_phosphorylation = 100.0*IL13_Rec*pJAK2;  // dwdp[8]
    dflux_v3_v_2_dRec_recycle = 100.0*Rec_i;  // dwdp[9]
    dflux_v14_v_13_dSOCS3_accumulation = -100.0*SOCS3_translation*SOCS3mRNA/std::pow(SOCS3_accumulation + SOCS3mRNA, 2);  // dwdp[10]
    dflux_v15_v_14_dSOCS3_degradation = 100.0*SOCS3;  // dwdp[11]
    dflux_v14_v_13_dSOCS3_translation = 100.0*SOCS3mRNA/(SOCS3_accumulation + SOCS3mRNA);  // dwdp[12]
    dflux_v13_v_12_dSOCS3mRNA_production = 10.0*pSTAT5;  // dwdp[13]
    dflux_v10_v_9_dSTAT5_phosphorylation = 100.0*STAT5*pJAK2;  // dwdp[14]
    dflux_v9_v_8_dpJAK2_dephosphorylation = 9100.0*pJAK2;  // dwdp[15]
    dflux_v8_v_7_dpRec_degradation = 100.0*p_IL13_Rec_i;  // dwdp[16]
    dflux_v7_v_6_dpRec_intern = 100.0*p_IL13_Rec;  // dwdp[17]
    dflux_v11_v_10_dpSTAT5_dephosphorylation = 9100.0*pSTAT5;  // dwdp[18]
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
