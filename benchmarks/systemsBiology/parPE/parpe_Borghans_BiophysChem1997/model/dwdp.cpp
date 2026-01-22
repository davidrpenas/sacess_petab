#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Borghans_BiophysChem1997 {

static constexpr std::array<sunindextype, 24> dwdp_colptrs_Borghans_BiophysChem1997_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 15, 15, 15, 16, 17, 18, 18, 18, 18
};

void dwdp_colptrs_Borghans_BiophysChem1997(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Borghans_BiophysChem1997_));
}
} // namespace model_Borghans_BiophysChem1997
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Borghans_BiophysChem1997 {

static constexpr std::array<sunindextype, 18> dwdp_rowvals_Borghans_BiophysChem1997_ = {
    1, 4, 2, 6, 3, 6, 2, 2, 6, 1, 2, 5, 0, 5, 7, 6, 0, 0
};

void dwdp_rowvals_Borghans_BiophysChem1997(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Borghans_BiophysChem1997_));
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
#include "dwdp.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void dwdp_Borghans_BiophysChem1997(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v6_v_5_dVp = 1.0*beta_par;  // dwdp[11]
        dflux_v1_v_0_dbeta_par = 1.0*v1;  // dwdp[12]
        dflux_v6_v_5_dbeta_par = 1.0*Vp;  // dwdp[13]
        dflux_v1_v_0_dv0 = 1.0;  // dwdp[16]
        dflux_v1_v_0_dv1 = 1.0*beta_par;  // dwdp[17]
    }

    // dynamic expressions
    dflux_v2_v_1_dK2 = -2.0*K2*Vm2*std::pow(Z_state, 2)/std::pow(std::pow(K2, 2) + std::pow(Z_state, 2), 2);  // dwdp[0]
    dflux_v5_v_4_dK_par = 1.0*Z_state;  // dwdp[1]
    dflux_v3_v_2_dKa = -4.0*std::pow(A_state, 4)*std::pow(Ka, 3)*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 4)/(std::pow(std::pow(A_state, 4) + std::pow(Ka, 4), 2)*(std::pow(Ky, 2) + std::pow(Y_state, 2))*(std::pow(Kz, 4) + std::pow(Z_state, 4)));  // dwdp[2]
    dflux_v7_v_6_dKd = -1.0*std::pow(A_state, 2)*std::pow(Kd, n_par - 1)*Vd*std::pow(Z_state, n_par)*n_par/((std::pow(A_state, 2) + std::pow(Kp, 2))*std::pow(std::pow(Kd, n_par) + std::pow(Z_state, n_par), 2));  // dwdp[3]
    dflux_v4_v_3_dKf = 1.0*Y_state;  // dwdp[4]
    dflux_v7_v_6_dKp = -2.0*std::pow(A_state, 2)*Kp*Vd*std::pow(Z_state, n_par)/(std::pow(std::pow(A_state, 2) + std::pow(Kp, 2), 2)*(std::pow(Kd, n_par) + std::pow(Z_state, n_par)));  // dwdp[5]
    dflux_v3_v_2_dKy = -2.0*std::pow(A_state, 4)*Ky*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 4)/((std::pow(A_state, 4) + std::pow(Ka, 4))*std::pow(std::pow(Ky, 2) + std::pow(Y_state, 2), 2)*(std::pow(Kz, 4) + std::pow(Z_state, 4)));  // dwdp[6]
    dflux_v3_v_2_dKz = -4.0*std::pow(A_state, 4)*std::pow(Kz, 3)*Vm3*std::pow(Y_state, 2)*std::pow(Z_state, 4)/((std::pow(A_state, 4) + std::pow(Ka, 4))*(std::pow(Ky, 2) + std::pow(Y_state, 2))*std::pow(std::pow(Kz, 4) + std::pow(Z_state, 4), 2));  // dwdp[7]
    dflux_v7_v_6_dVd = 1.0*std::pow(A_state, 2)*std::pow(Z_state, n_par)/((std::pow(A_state, 2) + std::pow(Kp, 2))*(std::pow(Kd, n_par) + std::pow(Z_state, n_par)));  // dwdp[8]
    dflux_v2_v_1_dVm2 = 1.0*std::pow(Z_state, 2)/(std::pow(K2, 2) + std::pow(Z_state, 2));  // dwdp[9]
    dflux_v3_v_2_dVm3 = 1.0*std::pow(A_state, 4)*std::pow(Y_state, 2)*std::pow(Z_state, 4)/((std::pow(A_state, 4) + std::pow(Ka, 4))*(std::pow(Ky, 2) + std::pow(Y_state, 2))*(std::pow(Kz, 4) + std::pow(Z_state, 4)));  // dwdp[10]
    dflux_v8_v_7_depsilon_par = 1.0*A_state;  // dwdp[14]
    dflux_v7_v_6_dn_par = -1.0*std::pow(A_state, 2)*Vd*std::pow(Z_state, n_par)*(((Kd == 0) ? (
       Kd
    )
    : (
       std::pow(Kd, n_par)*std::log(Kd)
    )) + ((Z_state == 0) ? (
       Z_state
    )
    : (
       std::pow(Z_state, n_par)*std::log(Z_state)
    )))/((std::pow(A_state, 2) + std::pow(Kp, 2))*std::pow(std::pow(Kd, n_par) + std::pow(Z_state, n_par), 2)) + 1.0*std::pow(A_state, 2)*Vd*((Z_state == 0) ? (
       Z_state
    )
    : (
       std::pow(Z_state, n_par)*std::log(Z_state)
    ))/((std::pow(A_state, 2) + std::pow(Kp, 2))*(std::pow(Kd, n_par) + std::pow(Z_state, n_par)));  // dwdp[15]
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
