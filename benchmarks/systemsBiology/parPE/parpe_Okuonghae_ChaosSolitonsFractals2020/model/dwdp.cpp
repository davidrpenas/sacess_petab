#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 17> dwdp_colptrs_Okuonghae_ChaosSolitonsFractals2020_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 13, 13, 13, 13, 13
};

void dwdp_colptrs_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Okuonghae_ChaosSolitonsFractals2020_));
}
} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

static constexpr std::array<sunindextype, 13> dwdp_rowvals_Okuonghae_ChaosSolitonsFractals2020_ = {
    1, 8, 7, 9, 5, 6, 10, 11, 3, 4, 3, 4, 1
};

void dwdp_rowvals_Okuonghae_ChaosSolitonsFractals2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Okuonghae_ChaosSolitonsFractals2020_));
}
} // namespace model_Okuonghae_ChaosSolitonsFractals2020
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
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void dwdp_Okuonghae_ChaosSolitonsFractals2020(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dlam_dtransmission_rate_effective = (1 - delta)*(1 - eps)*(alpha*asymptomatic + symptomatic)/(-detected + total_pop);  // dwdp[0]
    dflux_asymptomatic_to_recovered_dgamma_a = 1.0*asymptomatic;  // dwdp[1]
    dflux_symptomatic_to_recovered_dgamma_0 = 1.0*symptomatic;  // dwdp[2]
    dflux_detected_to_recovered_dgamma_i = 1.0*detected;  // dwdp[3]
    dflux_symptomatic_to_detected_dpsi = 1.0*symptomatic;  // dwdp[4]
    dflux_asymptomatic_to_detected_dtheta = 1.0*asymptomatic;  // dwdp[5]
    dflux_symptomatic_to_deceased_dd_0 = 1.0*symptomatic;  // dwdp[6]
    dflux_detected_to_deceased_dd_D = 1.0*detected;  // dwdp[7]
    dflux_exposed_to_symptomatic_dsigma = 1.0*exposed*(1 - nu);  // dwdp[8]
    dflux_exposed_to_asymptopmatic_dsigma = 1.0*exposed*nu;  // dwdp[9]
    dflux_exposed_to_symptomatic_dnu = -1.0*exposed*sigma;  // dwdp[10]
    dflux_exposed_to_asymptopmatic_dnu = 1.0*exposed*sigma;  // dwdp[11]
    dlam_dalpha = asymptomatic*transmission_rate_effective*(1 - delta)*(1 - eps)/(-detected + total_pop);  // dwdp[12]
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
