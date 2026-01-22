#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 15> dwdp_colptrs_Oliveira_NatCommun2021_ = {
    0, 1, 4, 7, 9, 10, 11, 12, 12, 12, 12, 12, 12, 12, 12
};

void dwdp_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Oliveira_NatCommun2021_));
}
} // namespace model_Oliveira_NatCommun2021
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 12> dwdp_rowvals_Oliveira_NatCommun2021_ = {
    10, 7, 8, 9, 3, 4, 5, 6, 19, 0, 0, 0
};

void dwdp_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Oliveira_NatCommun2021_));
}
} // namespace model_Oliveira_NatCommun2021
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
namespace model_Oliveira_NatCommun2021 {

void dwdp_Oliveira_NatCommun2021(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dsymptomatic_to_hospitalized_dh_hosp_rate = gamma_s*xi;  // dwdp[1]
        dsymptomatic_to_icu_dh_hosp_rate = gamma_s*(1 - xi);  // dwdp[2]
        dsymptomatic_to_recovered_dh_hosp_rate = -gamma_s;  // dwdp[3]
        dhospitalized_to_deaths_dgamma_h = mu_h*(1 - omega_h);  // dwdp[4]
        dhospitalized_to_icu_dgamma_h = omega_h;  // dwdp[5]
        dhospitalized_to_recovered_dgamma_h = (1 - mu_h)*(1 - omega_h);  // dwdp[6]
        dicu_to_deaths_dgamma_u = mu_u*(1 - omega_u);  // dwdp[7]
    }

    // dynamic expressions
    dflux_Susceptible_to_Exposed_ddelta_ = 1.0*Asymptomatic*Susceptible*beta/population;  // dwdp[0]
    dflux_ICU_to_Hospitalized_dgamma_u = 1.0*ICU*(mu_u*omega_u - mu_u + 1);  // dwdp[8]
    dbeta_dbeta_0 = 1 - Heaviside_0;  // dwdp[9]
    dbeta_dbeta_1 = Heaviside_0*(Heaviside_2*(1 - beta_2_multiplier) - Heaviside_2 + 1);  // dwdp[10]
    dbeta_dbeta_2 = Heaviside_0*Heaviside_2*beta_2_multiplier;  // dwdp[11]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
