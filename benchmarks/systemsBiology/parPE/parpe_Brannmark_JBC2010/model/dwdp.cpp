#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 22> dwdp_colptrs_Brannmark_JBC2010_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 14, 14, 14, 14, 14, 14
};

void dwdp_colptrs_Brannmark_JBC2010(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 14> dwdp_rowvals_Brannmark_JBC2010_ = {
    1, 1, 2, 3, 4, 5, 5, 6, 7, 8, 8, 10, 9, 11
};

void dwdp_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
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
namespace model_Brannmark_JBC2010 {

void dwdp_Brannmark_JBC2010(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_v1_v_0_dk1a = 1.0*IR*insulin;  // dwdp[0]
    dflux_v1_v_0_dk1aBasic = 1.0*IR;  // dwdp[1]
    dflux_v2_v_1_dk1b = 1.0*IRins;  // dwdp[2]
    dflux_v3_v_2_dk1c = 1.0*IRins;  // dwdp[3]
    dflux_v4_v_3_dk1d = 1.0*IRp;  // dwdp[4]
    dflux_v5_v_4_dk1e = 1.0*IRiP;  // dwdp[5]
    dflux_v5_v_4_dk1f = 1.0*IRiP*Xp/(Xp + 1);  // dwdp[6]
    dflux_v6_v_5_dk1g = 1.0*IRp;  // dwdp[7]
    dflux_v7_v_6_dk1r = 1.0*IRi;  // dwdp[8]
    dflux_v8_v_7_dk21 = 1.0*IRS*(IRiP*k22 + IRp);  // dwdp[9]
    dflux_v8_v_7_dk22 = 1.0*IRS*IRiP*k21;  // dwdp[10]
    dflux_v10_v_9_dk3 = 1.0*IRSiP*X;  // dwdp[11]
    dflux_v9_v_8_dkm2 = 1.0*IRSiP;  // dwdp[12]
    dflux_v11_v_10_dkm3 = 1.0*Xp;  // dwdp[13]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
