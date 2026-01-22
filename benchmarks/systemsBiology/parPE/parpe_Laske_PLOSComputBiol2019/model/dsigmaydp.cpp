#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void dsigmaydp_Laske_PLOSComputBiol2019(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 16:
            dsigmaydp[0] = 1;
            break;
        case 17:
            dsigmaydp[1] = 1;
            break;
        case 18:
            dsigmaydp[2] = 1;
            break;
        case 19:
            dsigmaydp[3] = 1;
            break;
        case 20:
            dsigmaydp[4] = 13.0334748500416;
            break;
        case 21:
            dsigmaydp[5] = 3.1095605262909198;
            break;
        case 22:
            dsigmaydp[6] = 3.4465091034262501;
            break;
        case 23:
            dsigmaydp[7] = 1.9654070825149701;
            break;
        case 24:
            dsigmaydp[8] = 3.0367677443843699;
            break;
        case 25:
            dsigmaydp[9] = 4.0031196168321896;
            break;
        case 26:
            dsigmaydp[10] = 4.8783911282307004;
            break;
        case 27:
            dsigmaydp[11] = 3.9635379818204601;
            break;
        case 28:
            dsigmaydp[12] = 4.8820897165046002;
            break;
    }
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
