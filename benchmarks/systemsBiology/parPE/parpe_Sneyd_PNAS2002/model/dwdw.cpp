#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 14> dwdw_colptrs_Sneyd_PNAS2002_ = {
    0, 4, 6, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9
};

void dwdw_colptrs_Sneyd_PNAS2002(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Sneyd_PNAS2002_));
}
} // namespace model_Sneyd_PNAS2002
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 9> dwdw_rowvals_Sneyd_PNAS2002_ = {
    4, 5, 8, 9, 4, 5, 3, 7, 11
};

void dwdw_rowvals_Sneyd_PNAS2002(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Sneyd_PNAS2002_));
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
#include "dwdw.h"

namespace amici {
namespace model_Sneyd_PNAS2002 {

void dwdw_Sneyd_PNAS2002(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_v_1_dL1 = Ca*IP3*IPR_R*L3*(Ca*l4 + L3*k2)/(std::pow(L1, 2)*std::pow(Ca*(1 + L3/L1) + L3, 2));  // dwdw[0]
    dflux_v_2_dL1 = Ca*IPR_R*k1/(Ca*(L1/L3 + 1) + L1) - Ca*IPR_R*(Ca/L3 + 1)*(L1*k1 + l2)/std::pow(Ca*(L1/L3 + 1) + L1, 2);  // dwdw[1]
    dflux_v_5_dL1 = -IPR_A*L1*(k_4 + l_6)/std::pow(Ca + L1, 2) + IPR_A*(k_4 + l_6)/(Ca + L1);  // dwdw[2]
    dflux_v_6_dL1 = Ca*IPR_A*k1/(Ca + L1) - Ca*IPR_A*(L1*k1 + l2)/std::pow(Ca + L1, 2);  // dwdw[3]
    dflux_v_1_dL3 = IP3*IPR_R*k2/(Ca*(1 + L3/L1) + L3) - IP3*IPR_R*(Ca/L1 + 1)*(Ca*l4 + L3*k2)/std::pow(Ca*(1 + L3/L1) + L3, 2);  // dwdw[4]
    dflux_v_2_dL3 = std::pow(Ca, 2)*IPR_R*L1*(L1*k1 + l2)/(std::pow(L3, 2)*std::pow(Ca*(L1/L3 + 1) + L1, 2));  // dwdw[5]
    dflux_v_0_dL5 = Ca*IPR_O*(Ca*l_4 + k_2)/(std::pow(L5, 2)*std::pow(Ca/L5 + 1, 2));  // dwdw[6]
    dflux_v_4_dL5 = Ca*IPR_O*k4/(Ca + L5) - Ca*IPR_O*(L5*k4 + l6)/std::pow(Ca + L5, 2);  // dwdw[7]
    dflux_v_8_dL5 = -IPR_O*L5*k3/std::pow(Ca + L5, 2) + IPR_O*k3/(Ca + L5);  // dwdw[8]
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
