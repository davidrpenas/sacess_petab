#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<sunindextype, 151> dwdp_colptrs_Fujita_SciSignal2010_ = {
    0, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14
};

void dwdp_colptrs_Fujita_SciSignal2010(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Fujita_SciSignal2010_));
}
} // namespace model_Fujita_SciSignal2010
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<sunindextype, 14> dwdp_rowvals_Fujita_SciSignal2010_ = {
    10, 11, 1, 1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 9
};

void dwdp_rowvals_Fujita_SciSignal2010(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Fujita_SciSignal2010_));
}
} // namespace model_Fujita_SciSignal2010
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
namespace model_Fujita_SciSignal2010 {

void dwdp_Fujita_SciSignal2010(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v11_reaction_11_dEGFR_turnover = 68190.0;  // dwdp[1]
    }

    // dynamic expressions
    dflux_v10_reaction_10_dEGFR_turnover = 1.0*EGFR;  // dwdp[0]
    dflux_v1_reaction_1_dreaction_1_k1 = 1.0*EGF*EGFR;  // dwdp[2]
    dflux_v1_reaction_1_dreaction_1_k2 = -1.0*EGF_EGFR;  // dwdp[3]
    dflux_v2_reaction_2_dreaction_2_k1 = 1.0*Akt*pEGFR;  // dwdp[4]
    dflux_v2_reaction_2_dreaction_2_k2 = -1.0*pEGFR_Akt;  // dwdp[5]
    dflux_v3_reaction_3_dreaction_3_k1 = 1.0*pEGFR_Akt;  // dwdp[6]
    dflux_v4_reaction_4_dreaction_4_k1 = 1.0*pEGFR;  // dwdp[7]
    dflux_v5_reaction_5_dreaction_5_k1 = 1.0*S6*pAkt;  // dwdp[8]
    dflux_v5_reaction_5_dreaction_5_k2 = -1.0*pAkt_S6;  // dwdp[9]
    dflux_v6_reaction_6_dreaction_6_k1 = 1.0*pAkt_S6;  // dwdp[10]
    dflux_v7_reaction_7_dreaction_7_k1 = 1.0*pAkt;  // dwdp[11]
    dflux_v8_reaction_8_dreaction_8_k1 = 1.0*pS6;  // dwdp[12]
    dflux_v9_reaction_9_dreaction_9_k1 = 1.0*EGF_EGFR;  // dwdp[13]
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
