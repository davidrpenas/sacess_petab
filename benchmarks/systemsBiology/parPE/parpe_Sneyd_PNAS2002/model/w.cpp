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
namespace model_Sneyd_PNAS2002 {

void w_Sneyd_PNAS2002(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        L1 = k_1*l2/(k1*l_2);  // w[0]
        L3 = k_2*l4/(k2*l_4);  // w[1]
        L5 = k_4*l6/(k4*l_6);  // w[2]
    }

    // dynamic expressions
    flux_v_0 = IPR_O*(Ca*l_4 + k_2)/(Ca/L5 + 1);  // w[3]
    flux_v_1 = IP3*IPR_R*(Ca*l4 + L3*k2)/(Ca*(1 + L3/L1) + L3);  // w[4]
    flux_v_2 = Ca*IPR_R*(L1*k1 + l2)/(Ca*(L1/L3 + 1) + L1);  // w[5]
    flux_v_3 = IPR_I1*(k_1 + l_2);  // w[6]
    flux_v_4 = Ca*IPR_O*(L5*k4 + l6)/(Ca + L5);  // w[7]
    flux_v_5 = IPR_A*L1*(k_4 + l_6)/(Ca + L1);  // w[8]
    flux_v_6 = Ca*IPR_A*(L1*k1 + l2)/(Ca + L1);  // w[9]
    flux_v_7 = IPR_I2*(k_1 + l_2);  // w[10]
    flux_v_8 = IPR_O*L5*k3/(Ca + L5);  // w[11]
    flux_v_9 = IPR_S*k_3;  // w[12]
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
