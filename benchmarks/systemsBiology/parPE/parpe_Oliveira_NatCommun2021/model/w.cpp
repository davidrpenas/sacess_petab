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

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void w_Oliveira_NatCommun2021(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        exposed_to_asymptomatic = kappa*(1 - p_symp_rate);  // w[1]
        exposed_to_symptomatic = kappa*p_symp_rate;  // w[2]
        hospitalized_to_deaths = gamma_h*mu_h*(1 - omega_h);  // w[3]
        hospitalized_to_icu = gamma_h*omega_h;  // w[4]
        hospitalized_to_recovered = gamma_h*(1 - mu_h)*(1 - omega_h);  // w[5]
        icu_to_deaths = gamma_u*mu_u*(1 - omega_u);  // w[6]
        symptomatic_to_hospitalized = gamma_s*h_hosp_rate*xi;  // w[7]
        symptomatic_to_icu = gamma_s*h_hosp_rate*(1 - xi);  // w[8]
        symptomatic_to_recovered = gamma_s*(1 - h_hosp_rate);  // w[9]
    }

    // dynamic expressions
    beta = Heaviside_0*(Heaviside_2*(beta_1*(1 - beta_2_multiplier) + beta_2*beta_2_multiplier) + beta_1*(1 - Heaviside_2)) + beta_0*(1 - Heaviside_0);  // w[0]
    flux_Susceptible_to_Exposed = 1.0*Susceptible*beta*(Asymptomatic*delta_ + Symptomatic)/population;  // w[10]
    flux_Exposed_to_Asymptomatic = 1.0*Exposed*exposed_to_asymptomatic;  // w[11]
    flux_Exposed_to_Symptomatic = 1.0*Exposed*exposed_to_symptomatic;  // w[12]
    flux_Asymptomatic_to_Recovered = 1.0*Asymptomatic*gamma_a;  // w[13]
    flux_Symptomatic_to_Recovered = 1.0*Symptomatic*symptomatic_to_recovered;  // w[14]
    flux_Symptomatic_to_Hospitalized = 1.0*Symptomatic*symptomatic_to_hospitalized;  // w[15]
    flux_Hospitalized_to_Recovered = 1.0*Hospital*hospitalized_to_recovered;  // w[16]
    flux_Symptomatic_to_ICU = 1.0*Symptomatic*symptomatic_to_icu;  // w[17]
    flux_Hospitalized_to_ICU = 1.0*Hospital*hospitalized_to_icu;  // w[18]
    flux_ICU_to_Hospitalized = 1.0*ICU*gamma_u*(mu_u*omega_u - mu_u + 1);  // w[19]
    flux_Hospitalized_to_Deaths = 1.0*Hospital*hospitalized_to_deaths;  // w[20]
    flux_ICU_to_Deaths = 1.0*ICU*icu_to_deaths;  // w[21]
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
