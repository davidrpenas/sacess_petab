#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void sigmay_Schwen_PONE2014(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_IR1 = noiseParameter1_observable_IR1;  // sigmay[0]
    sigma_observable_IR2 = noiseParameter1_observable_IR2;  // sigmay[1]
    sigma_observable_IRsum = noiseParameter1_observable_IRsum;  // sigmay[2]
    sigma_observable_Insulin = noiseParameter1_observable_Insulin;  // sigmay[3]
    sigma_prior_ini_R1 = noiseParameter1_prior_ini_R1;  // sigmay[4]
    sigma_prior_ini_R2fold = noiseParameter1_prior_ini_R2fold;  // sigmay[5]
    sigma_prior_ka1 = noiseParameter1_prior_ka1;  // sigmay[6]
    sigma_prior_ka2fold = noiseParameter1_prior_ka2fold;  // sigmay[7]
    sigma_prior_kd1 = noiseParameter1_prior_kd1;  // sigmay[8]
    sigma_prior_kd2fold = noiseParameter1_prior_kd2fold;  // sigmay[9]
}

} // namespace model_Schwen_PONE2014
} // namespace amici
