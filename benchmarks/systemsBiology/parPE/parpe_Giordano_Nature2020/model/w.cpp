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
namespace model_Giordano_Nature2020 {

void w_Giordano_Nature2020(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    CurrentDiagnosedInfected = Diagnosed + Recognized + Threatened;  // w[0]
    CurrentTotalInfected = Ailing + Diagnosed + Infected + Recognized + Threatened;  // w[1]
    alpha = Heaviside_0*alpha_0 + Heaviside_2*alpha_28*(1 - Heaviside_6) + Heaviside_4*alpha_4*(1 - Heaviside_0) + Heaviside_6*alpha_22*(1 - Heaviside_4) + alpha_50*(1 - Heaviside_2);  // w[2]
    beta = Heaviside_0*beta_0 + Heaviside_2*beta_22*(1 - Heaviside_4) + Heaviside_4*beta_4*(1 - Heaviside_0) + beta_50*(1 - Heaviside_2);  // w[3]
    delta = Heaviside_0*delta_0 + Heaviside_2*delta_22*(1 - Heaviside_4) + Heaviside_4*delta_4*(1 - Heaviside_0) + delta_50*(1 - Heaviside_2);  // w[4]
    epsilon = Heaviside_10*epsilon_12*(1 - Heaviside_8) + Heaviside_2*epsilon_38*(1 - Heaviside_10) + Heaviside_8*epsilon_0 + epsilon_50*(1 - Heaviside_2);  // w[5]
    eta = Heaviside_10*eta_22*(1 - Heaviside_4) + Heaviside_2*eta_38*(1 - Heaviside_10) + Heaviside_4*eta_0 + eta_50*(1 - Heaviside_2);  // w[6]
    gamma = Heaviside_0*gamma_0 + Heaviside_2*gamma_28*(1 - Heaviside_6) + Heaviside_4*gamma_4*(1 - Heaviside_0) + Heaviside_6*gamma_22*(1 - Heaviside_4) + gamma_50*(1 - Heaviside_2);  // w[7]
    kappa = Heaviside_10*kappa_22*(1 - Heaviside_4) + Heaviside_2*kappa_38*(1 - Heaviside_10) + Heaviside_4*kappa_0 + kappa_50*(1 - Heaviside_2);  // w[8]
    lam = Heaviside_2*lam_22*(1 - Heaviside_4) + Heaviside_4*lam_0 + lam_50*(1 - Heaviside_2);  // w[9]
    mu = Heaviside_2*mu_22*(1 - Heaviside_4) + Heaviside_4*mu_0 + mu_50*(1 - Heaviside_2);  // w[10]
    nu = Heaviside_2*nu_22*(1 - Heaviside_4) + Heaviside_4*nu_0 + nu_50*(1 - Heaviside_2);  // w[11]
    rho = Heaviside_10*rho_22*(1 - Heaviside_4) + Heaviside_2*rho_38*(1 - Heaviside_10) + Heaviside_4*rho_0 + rho_50*(1 - Heaviside_2);  // w[12]
    sigma = Heaviside_10*sigma_22*(1 - Heaviside_4) + Heaviside_2*sigma_38*(1 - Heaviside_10) + Heaviside_4*sigma_0 + sigma_50*(1 - Heaviside_2);  // w[13]
    xi = Heaviside_10*xi_22*(1 - Heaviside_4) + Heaviside_2*xi_38*(1 - Heaviside_10) + Heaviside_4*xi_0 + xi_50*(1 - Heaviside_2);  // w[14]
    zeta = Heaviside_10*zeta_22*(1 - Heaviside_4) + Heaviside_2*zeta_38*(1 - Heaviside_10) + Heaviside_4*zeta_0 + zeta_50*(1 - Heaviside_2);  // w[15]
    CumulativeInfected = CurrentTotalInfected + Extinct + Healed;  // w[16]
    flux_Susceptible_to_Infected = 1.0*Susceptible*(Ailing*gamma + Diagnosed*beta + Infected*alpha + Recognized*delta);  // w[17]
    flux_Infected_to_Diagnosed = 1.0*Infected*epsilon;  // w[18]
    flux_Infected_to_Ailing = 1.0*Infected*zeta;  // w[19]
    flux_Infected_to_Healed = 1.0*Infected*lam;  // w[20]
    flux_Diagnosed_to_Recognized = 1.0*Diagnosed*eta;  // w[21]
    flux_Diagnosed_to_Healed = 1.0*Diagnosed*rho;  // w[22]
    flux_Ailing_to_Recognised = 1.0*Ailing*theta;  // w[23]
    flux_Ailing_to_Healed = 1.0*Ailing*kappa;  // w[24]
    flux_Ailing_to_Threatened = 1.0*Ailing*mu;  // w[25]
    flux_Recognised_to_Threatened = 1.0*Recognized*nu;  // w[26]
    flux_Recognised_to_Healed = 1.0*Recognized*xi;  // w[27]
    flux_Threatened_to_Extinct = 1.0*Threatened*tau;  // w[28]
    flux_Threatened_to_Healed = 1.0*Threatened*sigma;  // w[29]
}

} // namespace model_Giordano_Nature2020
} // namespace amici
