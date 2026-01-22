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
namespace model_Rahman_MBS2016 {

void sigmay_Rahman_MBS2016(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_prevalence = noiseParameter1_observable_prevalence;  // sigmay[0]
}

} // namespace model_Rahman_MBS2016
} // namespace amici
