#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void dsigmaydp_Blasi_CellSystems2016(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 8:
            dsigmaydp[0] = 1;
            break;
        case 9:
            dsigmaydp[1] = 1;
            break;
        case 10:
            dsigmaydp[2] = 1;
            break;
        case 11:
            dsigmaydp[3] = 1;
            break;
        case 12:
            dsigmaydp[4] = 1;
            break;
        case 13:
            dsigmaydp[5] = 1;
            break;
        case 14:
            dsigmaydp[6] = 1;
            break;
        case 15:
            dsigmaydp[7] = 1;
            break;
        case 16:
            dsigmaydp[8] = 1;
            break;
        case 17:
            dsigmaydp[9] = 1;
            break;
        case 18:
            dsigmaydp[10] = 1;
            break;
        case 19:
            dsigmaydp[11] = 1;
            break;
        case 20:
            dsigmaydp[12] = 1;
            break;
        case 21:
            dsigmaydp[13] = 1;
            break;
        case 22:
            dsigmaydp[14] = 1;
            break;
        case 23:
            dsigmaydp[15] = 1;
            break;
    }
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
