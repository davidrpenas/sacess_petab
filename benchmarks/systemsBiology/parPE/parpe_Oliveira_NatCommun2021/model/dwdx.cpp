#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 10> dwdx_colptrs_Oliveira_NatCommun2021_ = {
    0, 2, 3, 5, 9, 9, 12, 14, 14, 14
};

void dwdx_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Oliveira_NatCommun2021_));
}
} // namespace model_Oliveira_NatCommun2021
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<sunindextype, 14> dwdx_rowvals_Oliveira_NatCommun2021_ = {
    11, 12, 10, 10, 13, 10, 14, 15, 17, 16, 18, 20, 19, 21
};

void dwdx_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Oliveira_NatCommun2021_));
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
#include "dwdx.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void dwdx_Oliveira_NatCommun2021(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_Exposed_to_Asymptomatic_dExposed = 1.0*exposed_to_asymptomatic;  // dwdx[0]
        dflux_Exposed_to_Symptomatic_dExposed = 1.0*exposed_to_symptomatic;  // dwdx[1]
        dflux_Asymptomatic_to_Recovered_dAsymptomatic = 1.0*gamma_a;  // dwdx[4]
        dflux_Symptomatic_to_Recovered_dSymptomatic = 1.0*symptomatic_to_recovered;  // dwdx[6]
        dflux_Symptomatic_to_Hospitalized_dSymptomatic = 1.0*symptomatic_to_hospitalized;  // dwdx[7]
        dflux_Symptomatic_to_ICU_dSymptomatic = 1.0*symptomatic_to_icu;  // dwdx[8]
        dflux_Hospitalized_to_Recovered_dHospital = 1.0*hospitalized_to_recovered;  // dwdx[9]
        dflux_Hospitalized_to_ICU_dHospital = 1.0*hospitalized_to_icu;  // dwdx[10]
        dflux_Hospitalized_to_Deaths_dHospital = 1.0*hospitalized_to_deaths;  // dwdx[11]
        dflux_ICU_to_Hospitalized_dICU = 1.0*gamma_u*(mu_u*omega_u - mu_u + 1);  // dwdx[12]
        dflux_ICU_to_Deaths_dICU = 1.0*icu_to_deaths;  // dwdx[13]
    }

    // dynamic expressions
    dflux_Susceptible_to_Exposed_dSusceptible = 1.0*beta*(Asymptomatic*delta_ + Symptomatic)/population;  // dwdx[2]
    dflux_Susceptible_to_Exposed_dAsymptomatic = 1.0*Susceptible*beta*delta_/population;  // dwdx[3]
    dflux_Susceptible_to_Exposed_dSymptomatic = 1.0*Susceptible*beta/population;  // dwdx[5]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
