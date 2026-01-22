#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<sunindextype, 10> dwdx_colptrs_Fujita_SciSignal2010_ = {
    0, 2, 4, 6, 7, 9, 10, 12, 13, 15
};

void dwdx_colptrs_Fujita_SciSignal2010(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Fujita_SciSignal2010_));
}
} // namespace model_Fujita_SciSignal2010
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fujita_SciSignal2010 {

static constexpr std::array<sunindextype, 15> dwdx_rowvals_Fujita_SciSignal2010_ = {
    1, 10, 2, 4, 2, 3, 2, 5, 7, 5, 5, 6, 8, 1, 9
};

void dwdx_rowvals_Fujita_SciSignal2010(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Fujita_SciSignal2010_));
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
#include "dwdx.h"

namespace amici {
namespace model_Fujita_SciSignal2010 {

void dwdx_Fujita_SciSignal2010(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v10_reaction_10_dEGFR = 1.0*EGFR_turnover;  // dwdx[1]
        dflux_v4_reaction_4_dpEGFR = 1.0*reaction_4_k1;  // dwdx[3]
        dflux_v2_reaction_2_dpEGFR_Akt = -1.0*reaction_2_k2;  // dwdx[4]
        dflux_v3_reaction_3_dpEGFR_Akt = 1.0*reaction_3_k1;  // dwdx[5]
        dflux_v7_reaction_7_dpAkt = 1.0*reaction_7_k1;  // dwdx[8]
        dflux_v5_reaction_5_dpAkt_S6 = -1.0*reaction_5_k2;  // dwdx[10]
        dflux_v6_reaction_6_dpAkt_S6 = 1.0*reaction_6_k1;  // dwdx[11]
        dflux_v8_reaction_8_dpS6 = 1.0*reaction_8_k1;  // dwdx[12]
        dflux_v1_reaction_1_dEGF_EGFR = -1.0*reaction_1_k2;  // dwdx[13]
        dflux_v9_reaction_9_dEGF_EGFR = 1.0*reaction_9_k1;  // dwdx[14]
    }

    // dynamic expressions
    dflux_v1_reaction_1_dEGFR = 1.0*EGF*reaction_1_k1;  // dwdx[0]
    dflux_v2_reaction_2_dpEGFR = 1.0*Akt*reaction_2_k1;  // dwdx[2]
    dflux_v2_reaction_2_dAkt = 1.0*pEGFR*reaction_2_k1;  // dwdx[6]
    dflux_v5_reaction_5_dpAkt = 1.0*S6*reaction_5_k1;  // dwdx[7]
    dflux_v5_reaction_5_dS6 = 1.0*pAkt*reaction_5_k1;  // dwdx[9]
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
