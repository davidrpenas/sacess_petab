#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void sigmay_Beer_MolBioSystems2014(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_Bacnorm = noiseParameter1_Bacnorm;  // sigmay[0]
    sigma_IndconcNormRange = noiseParameter1_IndconcNormRange;  // sigmay[1]
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
