#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void w_Borghans_BiophysChem1997(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        flux_v1_v_0 = 1.0*beta_par*v1 + 1.0*v0;  // w[0]
        flux_v6_v_5 = 1.0*Vp*beta_par;  // w[5]
    }

    // dynamic expressions
    flux_v2_v_1 = 1.0*Vm2*std::pow(Z_state, 2)/(std::pow(K2, 2) + std::pow(Z_state, 2));  // w[1]
    flux_v3_v_2 = 1.0*std::pow(A_state, 4)*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 4)/((std::pow(A_state, 4) + std::pow(Ka, 4))*(std::pow(Ky, 2) + std::pow(Y_state, 2))*(std::pow(Kz, 4) + std::pow(Z_state, 4)));  // w[2]
    flux_v4_v_3 = 1.0*Kf*Y_state;  // w[3]
    flux_v5_v_4 = 1.0*K_par*Z_state;  // w[4]
    flux_v7_v_6 = 1.0*std::pow(A_state, 2)*Vd*std::pow(Z_state, n_par)/((std::pow(A_state, 2) + std::pow(Kp, 2))*(std::pow(Kd, n_par) + std::pow(Z_state, n_par)));  // w[6]
    flux_v8_v_7 = 1.0*A_state*epsilon_par;  // w[7]
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
