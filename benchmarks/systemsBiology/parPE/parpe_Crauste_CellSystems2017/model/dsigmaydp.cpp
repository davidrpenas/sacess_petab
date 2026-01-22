#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void dsigmaydp_Crauste_CellSystems2017(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 12:
            dsigmaydp[0] = 1;
            break;
        case 13:
            dsigmaydp[1] = 1;
            break;
        case 14:
            dsigmaydp[2] = 1;
            break;
        case 15:
            dsigmaydp[3] = 1;
            break;
        case 16:
            dsigmaydp[4] = 1;
            break;
        case 17:
            dsigmaydp[5] = 1;
            break;
        case 18:
            dsigmaydp[6] = 1;
            break;
        case 19:
            dsigmaydp[7] = 1;
            break;
        case 20:
            dsigmaydp[8] = 1;
            break;
        case 21:
            dsigmaydp[9] = 1;
            break;
        case 22:
            dsigmaydp[10] = 1;
            break;
        case 23:
            dsigmaydp[11] = 1;
            break;
        case 24:
            dsigmaydp[12] = 1;
            break;
        case 25:
            dsigmaydp[13] = 1;
            break;
        case 26:
            dsigmaydp[14] = 1;
            break;
        case 27:
            dsigmaydp[15] = 1;
            break;
        case 28:
            dsigmaydp[16] = 1;
            break;
        case 29:
            dsigmaydp[17] = 1;
            break;
        case 30:
            dsigmaydp[18] = 1;
            break;
        case 31:
            dsigmaydp[19] = 1;
            break;
        case 32:
            dsigmaydp[20] = 1;
            break;
    }
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
