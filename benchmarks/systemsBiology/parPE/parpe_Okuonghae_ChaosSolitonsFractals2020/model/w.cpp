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
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void w_Okuonghae_ChaosSolitonsFractals2020(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    total_pop = asymptomatic + detected + exposed + recovered + susceptible + symptomatic;  // w[0]
    lam = transmission_rate_effective*(1 - delta)*(1 - eps)*(alpha*asymptomatic + symptomatic)/(-detected + total_pop);  // w[1]
    flux_susceptible_to_exposed = 1.0*lam*susceptible;  // w[2]
    flux_exposed_to_symptomatic = 1.0*exposed*sigma*(1 - nu);  // w[3]
    flux_exposed_to_asymptopmatic = 1.0*exposed*nu*sigma;  // w[4]
    flux_symptomatic_to_detected = 1.0*psi*symptomatic;  // w[5]
    flux_asymptomatic_to_detected = 1.0*asymptomatic*theta;  // w[6]
    flux_symptomatic_to_recovered = 1.0*gamma_0*symptomatic;  // w[7]
    flux_asymptomatic_to_recovered = 1.0*asymptomatic*gamma_a;  // w[8]
    flux_detected_to_recovered = 1.0*detected*gamma_i;  // w[9]
    flux_symptomatic_to_deceased = 1.0*d_0*symptomatic;  // w[10]
    flux_detected_to_deceased = 1.0*d_D*detected;  // w[11]
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
