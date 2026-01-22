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
namespace model_Fujita_SciSignal2010 {

void w_Fujita_SciSignal2010(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        flux_v11_reaction_11 = 68190.0*EGFR_turnover;  // w[11]
    }

    // dynamic expressions
    EGF = Heaviside_0*(EGF_0 + EGF_rate*t);  // w[0]
    flux_v1_reaction_1 = 1.0*EGF*EGFR*reaction_1_k1 - 1.0*EGF_EGFR*reaction_1_k2;  // w[1]
    flux_v2_reaction_2 = 1.0*Akt*pEGFR*reaction_2_k1 - 1.0*pEGFR_Akt*reaction_2_k2;  // w[2]
    flux_v3_reaction_3 = 1.0*pEGFR_Akt*reaction_3_k1;  // w[3]
    flux_v4_reaction_4 = 1.0*pEGFR*reaction_4_k1;  // w[4]
    flux_v5_reaction_5 = 1.0*S6*pAkt*reaction_5_k1 - 1.0*pAkt_S6*reaction_5_k2;  // w[5]
    flux_v6_reaction_6 = 1.0*pAkt_S6*reaction_6_k1;  // w[6]
    flux_v7_reaction_7 = 1.0*pAkt*reaction_7_k1;  // w[7]
    flux_v8_reaction_8 = 1.0*pS6*reaction_8_k1;  // w[8]
    flux_v9_reaction_9 = 1.0*EGF_EGFR*reaction_9_k1;  // w[9]
    flux_v10_reaction_10 = 1.0*EGFR*EGFR_turnover;  // w[10]
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
