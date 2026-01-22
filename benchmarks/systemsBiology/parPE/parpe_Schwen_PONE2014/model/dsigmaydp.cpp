#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dsigmaydp_Schwen_PONE2014(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 23:
            dsigmaydp[0] = 1;
            break;
        case 24:
            dsigmaydp[1] = 1;
            break;
        case 25:
            dsigmaydp[2] = 1;
            break;
        case 26:
            dsigmaydp[3] = 1;
            break;
        case 27:
            dsigmaydp[4] = 1;
            break;
        case 28:
            dsigmaydp[5] = 1;
            break;
        case 29:
            dsigmaydp[6] = 1;
            break;
        case 30:
            dsigmaydp[7] = 1;
            break;
        case 31:
            dsigmaydp[8] = 1;
            break;
        case 32:
            dsigmaydp[9] = 1;
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
