#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Weber_BMC2015 {

void dsigmaydp_Weber_BMC2015(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 31:
            dsigmaydp[0] = 1;
            break;
        case 32:
            dsigmaydp[1] = 1;
            break;
        case 33:
            dsigmaydp[2] = 1;
            break;
        case 34:
            dsigmaydp[3] = 1;
            break;
        case 35:
            dsigmaydp[4] = 1;
            break;
        case 36:
            dsigmaydp[5] = 1;
            break;
        case 37:
            dsigmaydp[6] = 1;
            break;
        case 38:
            dsigmaydp[7] = 1;
            break;
    }
}

} // namespace model_Weber_BMC2015
} // namespace amici
