#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void dsigmaydp_Zheng_PNAS2012(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 45:
            dsigmaydp[0] = 1;
            break;
        case 46:
            dsigmaydp[1] = 1;
            break;
        case 47:
            dsigmaydp[2] = 1;
            break;
        case 48:
            dsigmaydp[3] = 1;
            break;
        case 49:
            dsigmaydp[4] = 1;
            break;
        case 50:
            dsigmaydp[5] = 1;
            break;
        case 51:
            dsigmaydp[6] = 1;
            break;
        case 52:
            dsigmaydp[7] = 1;
            break;
        case 53:
            dsigmaydp[8] = 1;
            break;
        case 54:
            dsigmaydp[9] = 1;
            break;
        case 55:
            dsigmaydp[10] = 1;
            break;
        case 56:
            dsigmaydp[11] = 1;
            break;
        case 57:
            dsigmaydp[12] = 1;
            break;
        case 58:
            dsigmaydp[13] = 1;
            break;
        case 59:
            dsigmaydp[14] = 1;
            break;
    }
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
