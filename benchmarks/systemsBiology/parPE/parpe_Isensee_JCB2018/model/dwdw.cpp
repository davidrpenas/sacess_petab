#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 52> dwdw_colptrs_Isensee_JCB2018_ = {
    0, 2, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
};

void dwdw_colptrs_Isensee_JCB2018(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Isensee_JCB2018_));
}
} // namespace model_Isensee_JCB2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<sunindextype, 8> dwdw_rowvals_Isensee_JCB2018_ = {
    13, 15, 47, 12, 9, 10, 11, 8
};

void dwdw_rowvals_Isensee_JCB2018(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Isensee_JCB2018_));
}
} // namespace model_Isensee_JCB2018
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
#include "dwdw.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void dwdw_Isensee_JCB2018(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v5_v_4_dIBMXex = 1.0*ki_IBMX;  // dwdw[3]
        dflux_v2_v_1_dRp8_Br_cAMPS_pAB = 1.0*ki_Rp8_Br_cAMPS_pAB*xi_i_Rp8_Br_cAMPS_pAB;  // dwdw[4]
        dflux_v3_v_2_dRp8_pCPT_cAMPS_pAB = 1.0*ki_Rp8_pCPT_cAMPS_pAB*xi_i_Rp8_pCPT_cAMPS_pAB;  // dwdw[5]
        dflux_v4_v_3_dRp_cAMPS_pAB = 1.0*ki_Rp_cAMPS_pAB*xi_i_Rp_cAMPS_pAB;  // dwdw[6]
        dflux_v1_v_0_dSp8_Br_cAMPS_AM = 1.0*ki_Sp8_Br_cAMPS_AM*xi_i_Sp8_Br_cAMPS_AM;  // dwdw[7]
    }

    // dynamic expressions
    dflux_v6_v_5_dFsk = 1.0*AC*kf_Fsk;  // dwdw[0]
    dflux_v8_v_7_dFsk = 1.0*kf_Fsk*pAC;  // dwdw[1]
    dflux_v40_v_39_dH89 = 1.0*Csub*kf_H89;  // dwdw[2]
}

} // namespace model_Isensee_JCB2018
} // namespace amici
