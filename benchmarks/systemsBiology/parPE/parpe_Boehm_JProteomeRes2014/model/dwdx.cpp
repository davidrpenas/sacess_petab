#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<sunindextype, 9> dwdx_colptrs_Boehm_JProteomeRes2014_ = {
    0, 2, 4, 5, 6, 7, 8, 9, 10
};

void dwdx_colptrs_Boehm_JProteomeRes2014(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Boehm_JProteomeRes2014_));
}
} // namespace model_Boehm_JProteomeRes2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<sunindextype, 10> dwdx_rowvals_Boehm_JProteomeRes2014_ = {
    1, 2, 2, 3, 5, 4, 6, 7, 8, 9
};

void dwdx_rowvals_Boehm_JProteomeRes2014(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Boehm_JProteomeRes2014_));
}
} // namespace model_Boehm_JProteomeRes2014
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
namespace model_Boehm_JProteomeRes2014 {

void dwdx_Boehm_JProteomeRes2014(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v5_v_4_dpApB = 1.3999999999999999*k_imp_hetero;  // dwdx[4]
        dflux_v4_v_3_dpApA = 1.3999999999999999*k_imp_homo;  // dwdx[5]
        dflux_v6_v_5_dpBpB = 1.3999999999999999*k_imp_homo;  // dwdx[6]
        dflux_v7_v_6_dnucpApA = 0.45000000000000001*k_exp_homo;  // dwdx[7]
        dflux_v8_v_7_dnucpApB = 0.45000000000000001*k_exp_hetero;  // dwdx[8]
        dflux_v9_v_8_dnucpBpB = 0.45000000000000001*k_exp_homo;  // dwdx[9]
    }

    // dynamic expressions
    dflux_v1_v_0_dSTAT5A = 2.7999999999999998*BaF3_Epo*STAT5A*k_phos;  // dwdx[0]
    dflux_v2_v_1_dSTAT5A = 1.3999999999999999*BaF3_Epo*STAT5B*k_phos;  // dwdx[1]
    dflux_v2_v_1_dSTAT5B = 1.3999999999999999*BaF3_Epo*STAT5A*k_phos;  // dwdx[2]
    dflux_v3_v_2_dSTAT5B = 2.7999999999999998*BaF3_Epo*STAT5B*k_phos;  // dwdx[3]
}

} // namespace model_Boehm_JProteomeRes2014
} // namespace amici
