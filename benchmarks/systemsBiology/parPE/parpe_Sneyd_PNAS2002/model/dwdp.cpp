#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 16> dwdp_colptrs_Sneyd_PNAS2002_ = {
    0, 3, 5, 6, 8, 11, 13, 14, 16, 19, 21, 23, 26, 28, 30, 30
};

void dwdp_colptrs_Sneyd_PNAS2002(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Sneyd_PNAS2002_));
}
} // namespace model_Sneyd_PNAS2002
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Sneyd_PNAS2002 {

static constexpr std::array<sunindextype, 30> dwdp_rowvals_Sneyd_PNAS2002_ = {
    0, 5, 9, 1, 4, 11, 2, 7, 0, 6, 10, 1, 3, 12, 2, 8, 0, 5, 9, 1, 4, 2, 7, 0, 6, 10, 1, 3, 2, 8
};

void dwdp_rowvals_Sneyd_PNAS2002(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Sneyd_PNAS2002_));
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
#include "dwdp.h"

namespace amici {
namespace model_Sneyd_PNAS2002 {

void dwdp_Sneyd_PNAS2002(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dL1_dk1 = -k_1*l2/(std::pow(k1, 2)*l_2);  // dwdp[0]
        dL3_dk2 = -k_2*l4/(std::pow(k2, 2)*l_4);  // dwdp[3]
        dL5_dk4 = -k_4*l6/(std::pow(k4, 2)*l_6);  // dwdp[6]
        dL1_dk_1 = l2/(k1*l_2);  // dwdp[8]
        dL3_dk_2 = l4/(k2*l_4);  // dwdp[11]
        dL5_dk_4 = l6/(k4*l_6);  // dwdp[14]
        dL1_dl2 = k_1/(k1*l_2);  // dwdp[16]
        dL3_dl4 = k_2/(k2*l_4);  // dwdp[19]
        dL5_dl6 = k_4/(k4*l_6);  // dwdp[21]
        dL1_dl_2 = -k_1*l2/(k1*std::pow(l_2, 2));  // dwdp[23]
        dL3_dl_4 = -k_2*l4/(k2*std::pow(l_4, 2));  // dwdp[26]
        dL5_dl_6 = -k_4*l6/(k4*std::pow(l_6, 2));  // dwdp[28]
    }

    // dynamic expressions
    dflux_v_2_dk1 = Ca*IPR_R*L1/(Ca*(L1/L3 + 1) + L1);  // dwdp[1]
    dflux_v_6_dk1 = Ca*IPR_A*L1/(Ca + L1);  // dwdp[2]
    dflux_v_1_dk2 = IP3*IPR_R*L3/(Ca*(1 + L3/L1) + L3);  // dwdp[4]
    dflux_v_8_dk3 = IPR_O*L5/(Ca + L5);  // dwdp[5]
    dflux_v_4_dk4 = Ca*IPR_O*L5/(Ca + L5);  // dwdp[7]
    dflux_v_3_dk_1 = IPR_I1;  // dwdp[9]
    dflux_v_7_dk_1 = IPR_I2;  // dwdp[10]
    dflux_v_0_dk_2 = IPR_O/(Ca/L5 + 1);  // dwdp[12]
    dflux_v_9_dk_3 = IPR_S;  // dwdp[13]
    dflux_v_5_dk_4 = IPR_A*L1/(Ca + L1);  // dwdp[15]
    dflux_v_2_dl2 = Ca*IPR_R/(Ca*(L1/L3 + 1) + L1);  // dwdp[17]
    dflux_v_6_dl2 = Ca*IPR_A/(Ca + L1);  // dwdp[18]
    dflux_v_1_dl4 = Ca*IP3*IPR_R/(Ca*(1 + L3/L1) + L3);  // dwdp[20]
    dflux_v_4_dl6 = Ca*IPR_O/(Ca + L5);  // dwdp[22]
    dflux_v_3_dl_2 = IPR_I1;  // dwdp[24]
    dflux_v_7_dl_2 = IPR_I2;  // dwdp[25]
    dflux_v_0_dl_4 = Ca*IPR_O/(Ca/L5 + 1);  // dwdp[27]
    dflux_v_5_dl_6 = IPR_A*L1/(Ca + L1);  // dwdp[29]
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
