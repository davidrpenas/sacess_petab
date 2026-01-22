#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void w_Raia_CancerResearch2011(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    flux_v1_v_0 = 226.5*Kon_IL13Rec*Rec*il13_level;  // w[0]
    flux_v2_v_1 = 100.0*Rec*Rec_intern;  // w[1]
    flux_v3_v_2 = 100.0*Rec_i*Rec_recycle;  // w[2]
    flux_v4_v_3 = 100.0*IL13_Rec*Rec_phosphorylation*pJAK2;  // w[3]
    flux_v5_v_4 = 100.0*IL13_Rec*JAK2*JAK2_phosphorylation/(JAK2_p_inhibition*SOCS3 + 1);  // w[4]
    flux_v6_v_5 = 100.0*JAK2*JAK2_phosphorylation*p_IL13_Rec/(JAK2_p_inhibition*SOCS3 + 1);  // w[5]
    flux_v7_v_6 = 100.0*pRec_intern*p_IL13_Rec;  // w[6]
    flux_v8_v_7 = 100.0*pRec_degradation*p_IL13_Rec_i;  // w[7]
    flux_v9_v_8 = 9100.0*pJAK2*pJAK2_dephosphorylation;  // w[8]
    flux_v10_v_9 = 100.0*STAT5*STAT5_phosphorylation*pJAK2;  // w[9]
    flux_v11_v_10 = 9100.0*pSTAT5*pSTAT5_dephosphorylation;  // w[10]
    flux_v12_v_11 = 226.5*DecoyR*DecoyR_binding*il13_level;  // w[11]
    flux_v13_v_12 = 10.0*SOCS3mRNA_production*pSTAT5;  // w[12]
    flux_v14_v_13 = 100.0*SOCS3_translation*SOCS3mRNA/(SOCS3_accumulation + SOCS3mRNA);  // w[13]
    flux_v15_v_14 = 100.0*SOCS3*SOCS3_degradation;  // w[14]
    flux_v16_v_15 = 10.0*CD274mRNA_production*pSTAT5;  // w[15]
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
