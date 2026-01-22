#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<sunindextype, 15> dwdx_colptrs_Raia_CancerResearch2011_ = {
    0, 2, 3, 5, 7, 8, 10, 13, 14, 17, 18, 19, 19, 22, 22
};

void dwdx_colptrs_Raia_CancerResearch2011(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Raia_CancerResearch2011_));
}
} // namespace model_Raia_CancerResearch2011
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<sunindextype, 22> dwdx_rowvals_Raia_CancerResearch2011_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 4, 5, 3, 8, 9, 9, 10, 12, 15, 13, 11, 4, 5, 14
};

void dwdx_rowvals_Raia_CancerResearch2011(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Raia_CancerResearch2011_));
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
#include "dwdx.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void dwdx_Raia_CancerResearch2011(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v1_v_0_dRec = 226.5*Kon_IL13Rec*il13_level;  // dwdx[0]
        dflux_v2_v_1_dRec = 100.0*Rec_intern;  // dwdx[1]
        dflux_v3_v_2_dRec_i = 100.0*Rec_recycle;  // dwdx[2]
        dflux_v7_v_6_dp_IL13_Rec = 100.0*pRec_intern;  // dwdx[6]
        dflux_v8_v_7_dp_IL13_Rec_i = 100.0*pRec_degradation;  // dwdx[7]
        dflux_v9_v_8_dpJAK2 = 9100.0*pJAK2_dephosphorylation;  // dwdx[11]
        dflux_v11_v_10_dpSTAT5 = 9100.0*pSTAT5_dephosphorylation;  // dwdx[14]
        dflux_v13_v_12_dpSTAT5 = 10.0*SOCS3mRNA_production;  // dwdx[15]
        dflux_v16_v_15_dpSTAT5 = 10.0*CD274mRNA_production;  // dwdx[16]
        dflux_v12_v_11_dDecoyR = 226.5*DecoyR_binding*il13_level;  // dwdx[18]
        dflux_v15_v_14_dSOCS3 = 100.0*SOCS3_degradation;  // dwdx[21]
    }

    // dynamic expressions
    dflux_v4_v_3_dIL13_Rec = 100.0*Rec_phosphorylation*pJAK2;  // dwdx[3]
    dflux_v5_v_4_dIL13_Rec = 100.0*JAK2*JAK2_phosphorylation/(JAK2_p_inhibition*SOCS3 + 1);  // dwdx[4]
    dflux_v6_v_5_dp_IL13_Rec = 100.0*JAK2*JAK2_phosphorylation/(JAK2_p_inhibition*SOCS3 + 1);  // dwdx[5]
    dflux_v5_v_4_dJAK2 = 100.0*IL13_Rec*JAK2_phosphorylation/(JAK2_p_inhibition*SOCS3 + 1);  // dwdx[8]
    dflux_v6_v_5_dJAK2 = 100.0*JAK2_phosphorylation*p_IL13_Rec/(JAK2_p_inhibition*SOCS3 + 1);  // dwdx[9]
    dflux_v4_v_3_dpJAK2 = 100.0*IL13_Rec*Rec_phosphorylation;  // dwdx[10]
    dflux_v10_v_9_dpJAK2 = 100.0*STAT5*STAT5_phosphorylation;  // dwdx[12]
    dflux_v10_v_9_dSTAT5 = 100.0*STAT5_phosphorylation*pJAK2;  // dwdx[13]
    dflux_v14_v_13_dSOCS3mRNA = -100.0*SOCS3_translation*SOCS3mRNA/std::pow(SOCS3_accumulation + SOCS3mRNA, 2) + 100.0*SOCS3_translation/(SOCS3_accumulation + SOCS3mRNA);  // dwdx[17]
    dflux_v5_v_4_dSOCS3 = -100.0*IL13_Rec*JAK2*JAK2_p_inhibition*JAK2_phosphorylation/std::pow(JAK2_p_inhibition*SOCS3 + 1, 2);  // dwdx[19]
    dflux_v6_v_5_dSOCS3 = -100.0*JAK2*JAK2_p_inhibition*JAK2_phosphorylation*p_IL13_Rec/std::pow(JAK2_p_inhibition*SOCS3 + 1, 2);  // dwdx[20]
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
