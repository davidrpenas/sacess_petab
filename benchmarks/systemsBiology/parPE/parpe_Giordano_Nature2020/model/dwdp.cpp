#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 51> dwdp_colptrs_Giordano_Nature2020_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 43, 43, 43, 43, 43, 43, 43
};

void dwdp_colptrs_Giordano_Nature2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Giordano_Nature2020_));
}
} // namespace model_Giordano_Nature2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 43> dwdp_rowvals_Giordano_Nature2020_ = {
    23, 28, 2, 2, 2, 3, 3, 7, 7, 7, 4, 4, 5, 5, 10, 11, 12, 12, 8, 8, 14, 14, 13, 13, 15, 15, 6, 6, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 9, 13, 14, 15, 9
};

void dwdp_rowvals_Giordano_Nature2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Giordano_Nature2020_));
}
} // namespace model_Giordano_Nature2020
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
namespace model_Giordano_Nature2020 {

void dwdp_Giordano_Nature2020(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_Ailing_to_Recognised_dtheta = 1.0*Ailing;  // dwdp[0]
    dflux_Threatened_to_Extinct_dtau = 1.0*Threatened;  // dwdp[1]
    dalpha_dalpha_4 = Heaviside_4*(1 - Heaviside_0);  // dwdp[2]
    dalpha_dalpha_22 = Heaviside_6*(1 - Heaviside_4);  // dwdp[3]
    dalpha_dalpha_28 = Heaviside_2*(1 - Heaviside_6);  // dwdp[4]
    dbeta_dbeta_4 = Heaviside_4*(1 - Heaviside_0);  // dwdp[5]
    dbeta_dbeta_22 = Heaviside_2*(1 - Heaviside_4);  // dwdp[6]
    dgamma_dgamma_4 = Heaviside_4*(1 - Heaviside_0);  // dwdp[7]
    dgamma_dgamma_22 = Heaviside_6*(1 - Heaviside_4);  // dwdp[8]
    dgamma_dgamma_28 = Heaviside_2*(1 - Heaviside_6);  // dwdp[9]
    ddelta_ddelta_4 = Heaviside_4*(1 - Heaviside_0);  // dwdp[10]
    ddelta_ddelta_22 = Heaviside_2*(1 - Heaviside_4);  // dwdp[11]
    depsilon_depsilon_12 = Heaviside_10*(1 - Heaviside_8);  // dwdp[12]
    depsilon_depsilon_38 = Heaviside_2*(1 - Heaviside_10);  // dwdp[13]
    dmu_dmu_22 = Heaviside_2*(1 - Heaviside_4);  // dwdp[14]
    dnu_dnu_22 = Heaviside_2*(1 - Heaviside_4);  // dwdp[15]
    drho_drho_22 = Heaviside_10*(1 - Heaviside_4);  // dwdp[16]
    drho_drho_38 = Heaviside_2*(1 - Heaviside_10);  // dwdp[17]
    dkappa_dkappa_22 = Heaviside_10*(1 - Heaviside_4);  // dwdp[18]
    dkappa_dkappa_38 = Heaviside_2*(1 - Heaviside_10);  // dwdp[19]
    dxi_dxi_22 = Heaviside_10*(1 - Heaviside_4);  // dwdp[20]
    dxi_dxi_38 = Heaviside_2*(1 - Heaviside_10);  // dwdp[21]
    dsigma_dsigma_22 = Heaviside_10*(1 - Heaviside_4);  // dwdp[22]
    dsigma_dsigma_38 = Heaviside_2*(1 - Heaviside_10);  // dwdp[23]
    dzeta_dzeta_22 = Heaviside_10*(1 - Heaviside_4);  // dwdp[24]
    dzeta_dzeta_38 = Heaviside_2*(1 - Heaviside_10);  // dwdp[25]
    deta_deta_22 = Heaviside_10*(1 - Heaviside_4);  // dwdp[26]
    deta_deta_38 = Heaviside_2*(1 - Heaviside_10);  // dwdp[27]
    dalpha_dalpha_0 = Heaviside_0;  // dwdp[28]
    dbeta_dbeta_0 = Heaviside_0;  // dwdp[29]
    ddelta_ddelta_0 = Heaviside_0;  // dwdp[30]
    depsilon_depsilon_0 = Heaviside_8;  // dwdp[31]
    deta_deta_0 = Heaviside_4;  // dwdp[32]
    dgamma_dgamma_0 = Heaviside_0;  // dwdp[33]
    dkappa_dkappa_0 = Heaviside_4;  // dwdp[34]
    dmu_dmu_0 = Heaviside_4;  // dwdp[35]
    dnu_dnu_0 = Heaviside_4;  // dwdp[36]
    drho_drho_0 = Heaviside_4;  // dwdp[37]
    dlam_dlam_0 = Heaviside_4;  // dwdp[38]
    dsigma_dsigma_0 = Heaviside_4;  // dwdp[39]
    dxi_dxi_0 = Heaviside_4;  // dwdp[40]
    dzeta_dzeta_0 = Heaviside_4;  // dwdp[41]
    dlam_dlam_22 = Heaviside_2*(1 - Heaviside_4);  // dwdp[42]
}

} // namespace model_Giordano_Nature2020
} // namespace amici
