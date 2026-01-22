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
namespace model_Weber_BMC2015 {

void sigmay_Weber_BMC2015(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_yCERTpRN24 = noiseParameter1_yCERTpRN24;  // sigmay[0]
    sigma_yCERTt = noiseParameter1_yCERTt;  // sigmay[1]
    sigma_yPI4K3BpRN24 = noiseParameter1_yPI4K3BpRN24;  // sigmay[2]
    sigma_yPI4K3Bt = noiseParameter1_yPI4K3Bt;  // sigmay[3]
    sigma_yPKDpN0 = noiseParameter1_yPKDpN0;  // sigmay[4]
    sigma_yPKDpN24 = noiseParameter1_yPKDpN24;  // sigmay[5]
    sigma_yPKDpN25 = noiseParameter1_yPKDpN25;  // sigmay[6]
    sigma_yPKDt = noiseParameter1_yPKDt;  // sigmay[7]
}

} // namespace model_Weber_BMC2015
} // namespace amici
