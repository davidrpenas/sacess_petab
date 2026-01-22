#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void xdot_Okuonghae_ChaosSolitonsFractals2020(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dsusceptibledt = -flux_susceptible_to_exposed;  // xdot[0]
    dexposeddt = -flux_exposed_to_asymptopmatic - flux_exposed_to_symptomatic + flux_susceptible_to_exposed;  // xdot[1]
    dasymptomaticdt = -flux_asymptomatic_to_detected - flux_asymptomatic_to_recovered + flux_exposed_to_asymptopmatic;  // xdot[2]
    dsymptomaticdt = flux_exposed_to_symptomatic - flux_symptomatic_to_deceased - flux_symptomatic_to_detected - flux_symptomatic_to_recovered;  // xdot[3]
    ddetecteddt = flux_asymptomatic_to_detected - flux_detected_to_deceased - flux_detected_to_recovered + flux_symptomatic_to_detected;  // xdot[4]
    drecovereddt = flux_asymptomatic_to_recovered + flux_detected_to_recovered + flux_symptomatic_to_recovered;  // xdot[5]
    ddeceaseddt = flux_detected_to_deceased + flux_symptomatic_to_deceased;  // xdot[6]
    ddetected_cumulativedt = flux_asymptomatic_to_detected + flux_symptomatic_to_detected;  // xdot[7]
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
