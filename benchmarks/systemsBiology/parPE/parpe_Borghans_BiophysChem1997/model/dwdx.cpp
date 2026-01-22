#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Borghans_BiophysChem1997 {

static constexpr std::array<sunindextype, 4> dwdx_colptrs_Borghans_BiophysChem1997_ = {
    0, 4, 6, 9
};

void dwdx_colptrs_Borghans_BiophysChem1997(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Borghans_BiophysChem1997_));
}
} // namespace model_Borghans_BiophysChem1997
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Borghans_BiophysChem1997 {

static constexpr std::array<sunindextype, 9> dwdx_rowvals_Borghans_BiophysChem1997_ = {
    1, 2, 4, 6, 2, 3, 2, 6, 7
};

void dwdx_rowvals_Borghans_BiophysChem1997(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Borghans_BiophysChem1997_));
}
} // namespace model_Borghans_BiophysChem1997
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void dwdx_Borghans_BiophysChem1997(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v5_v_4_dZ_state = 1.0*K_par;  // dwdx[2]
        dflux_v4_v_3_dY_state = 1.0*Kf;  // dwdx[5]
        dflux_v8_v_7_dA_state = 1.0*epsilon_par;  // dwdx[8]
    }

    // dynamic expressions
    dflux_v2_v_1_dZ_state = -2.0*Vm2*std::pow(Z_state, 3)/std::pow(std::pow(K2, 2) + std::pow(Z_state, 2), 2) + 2.0*Vm2*Z_state/(std::pow(K2, 2) + std::pow(Z_state, 2));  // dwdx[0]
    dflux_v3_v_2_dZ_state = -4.0*std::pow(A_state, 4)*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 7)/((std::pow(A_state, 4) + std::pow(Ka, 4))*(std::pow(Ky, 2) + std::pow(Y_state, 2))*std::pow(std::pow(Kz, 4) + std::pow(Z_state, 4), 2)) + 4.0*std::pow(A_state, 4)*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 3)/((std::pow(A_state, 4) + std::pow(Ka, 4))*(std::pow(Ky, 2) + std::pow(Y_state, 2))*(std::pow(Kz, 4) + std::pow(Z_state, 4)));  // dwdx[1]
    dflux_v7_v_6_dZ_state = 1.0*std::pow(A_state, 2)*Vd*std::pow(Z_state, n_par - 1)*n_par/((std::pow(A_state, 2) + std::pow(Kp, 2))*(std::pow(Kd, n_par) + std::pow(Z_state, n_par))) - 1.0*std::pow(A_state, 2)*Vd*std::pow(Z_state, 2*n_par - 1)*n_par/((std::pow(A_state, 2) + std::pow(Kp, 2))*std::pow(std::pow(Kd, n_par) + std::pow(Z_state, n_par), 2));  // dwdx[3]
    dflux_v3_v_2_dY_state = -2.0*std::pow(A_state, 4)*Vm3*std::pow(Y_state, 3)*std::pow(Z_state, 4)/((std::pow(A_state, 4) + std::pow(Ka, 4))*std::pow(std::pow(Ky, 2) + std::pow(Y_state, 2), 2)*(std::pow(Kz, 4) + std::pow(Z_state, 4))) + 2.0*std::pow(A_state, 4)*Vm3*Y_state*std::pow(Z_state, 4)/((std::pow(A_state, 4) + std::pow(Ka, 4))*(std::pow(Ky, 2) + std::pow(Y_state, 2))*(std::pow(Kz, 4) + std::pow(Z_state, 4)));  // dwdx[4]
    dflux_v3_v_2_dA_state = -4.0*std::pow(A_state, 7)*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 4)/(std::pow(std::pow(A_state, 4) + std::pow(Ka, 4), 2)*(std::pow(Ky, 2) + std::pow(Y_state, 2))*(std::pow(Kz, 4) + std::pow(Z_state, 4))) + 4.0*std::pow(A_state, 3)*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 4)/((std::pow(A_state, 4) + std::pow(Ka, 4))*(std::pow(Ky, 2) + std::pow(Y_state, 2))*(std::pow(Kz, 4) + std::pow(Z_state, 4)));  // dwdx[6]
    dflux_v7_v_6_dA_state = -2.0*std::pow(A_state, 3)*Vd*std::pow(Z_state, n_par)/(std::pow(std::pow(A_state, 2) + std::pow(Kp, 2), 2)*(std::pow(Kd, n_par) + std::pow(Z_state, n_par))) + 2.0*A_state*Vd*std::pow(Z_state, n_par)/((std::pow(A_state, 2) + std::pow(Kp, 2))*(std::pow(Kd, n_par) + std::pow(Z_state, n_par)));  // dwdx[7]
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
