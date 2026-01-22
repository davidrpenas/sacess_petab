#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 7> dwdx_colptrs_Sneyd_PNAS2002_ = {
    0, 3, 5, 6, 7, 9, 10
};

void dwdx_colptrs_Sneyd_PNAS2002(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Sneyd_PNAS2002_));
}
} // namespace model_Sneyd_PNAS2002
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 10> dwdx_rowvals_Sneyd_PNAS2002_ = {
    3, 7, 11, 4, 5, 6, 12, 8, 9, 10
};

void dwdx_rowvals_Sneyd_PNAS2002(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Sneyd_PNAS2002_));
}
} // namespace model_Sneyd_PNAS2002
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
namespace model_Sneyd_PNAS2002 {

void dwdx_Sneyd_PNAS2002(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v_0_dIPR_O = (Ca*l_4 + k_2)/(Ca/L5 + 1);  // dwdx[0]
        dflux_v_4_dIPR_O = Ca*(L5*k4 + l6)/(Ca + L5);  // dwdx[1]
        dflux_v_8_dIPR_O = L5*k3/(Ca + L5);  // dwdx[2]
        dflux_v_1_dIPR_R = IP3*(Ca*l4 + L3*k2)/(Ca*(1 + L3/L1) + L3);  // dwdx[3]
        dflux_v_2_dIPR_R = Ca*(L1*k1 + l2)/(Ca*(L1/L3 + 1) + L1);  // dwdx[4]
        dflux_v_3_dIPR_I1 = k_1 + l_2;  // dwdx[5]
        dflux_v_9_dIPR_S = k_3;  // dwdx[6]
        dflux_v_5_dIPR_A = L1*(k_4 + l_6)/(Ca + L1);  // dwdx[7]
        dflux_v_6_dIPR_A = Ca*(L1*k1 + l2)/(Ca + L1);  // dwdx[8]
        dflux_v_7_dIPR_I2 = k_1 + l_2;  // dwdx[9]
    }
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
