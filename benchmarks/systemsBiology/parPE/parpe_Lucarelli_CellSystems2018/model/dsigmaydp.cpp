#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void dsigmaydp_Lucarelli_CellSystems2018(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 72:
            dsigmaydp[0] = 1;
            break;
        case 73:
            dsigmaydp[1] = 1;
            break;
        case 74:
            dsigmaydp[2] = 1;
            break;
        case 75:
            dsigmaydp[3] = 1;
            break;
        case 76:
            dsigmaydp[4] = 1;
            break;
        case 77:
            dsigmaydp[5] = 1;
            break;
        case 78:
            dsigmaydp[6] = 1;
            break;
        case 79:
            dsigmaydp[59] = 1;
            break;
        case 80:
            dsigmaydp[60] = 1;
            break;
        case 81:
            dsigmaydp[61] = 1;
            break;
        case 82:
            dsigmaydp[62] = 1;
            break;
        case 83:
            dsigmaydp[64] = 1;
            break;
    }
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
