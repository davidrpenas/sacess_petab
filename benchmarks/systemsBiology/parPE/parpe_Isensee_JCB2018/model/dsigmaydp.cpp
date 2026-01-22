#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void dsigmaydp_Isensee_JCB2018(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 35:
            dsigmaydp[0] = noiseParameter1_pRII_Microscopy*rho_pRII_Microscopy;
            break;
        case 37:
            dsigmaydp[1] = rho_pRII_Western;
            break;
        case 40:
            dsigmaydp[2] = rho_Calpha_Microscopy;
            break;
        case 45:
            dsigmaydp[0] = rho_pRII_Microscopy*s_pRII_global;
            break;
        case 46:
            dsigmaydp[0] = noiseParameter1_pRII_Microscopy*s_pRII_global;
            break;
        case 47:
            dsigmaydp[1] = s_pRII_Western;
            break;
        case 48:
            dsigmaydp[2] = s_Calpha_global;
            break;
        case 49:
            dsigmaydp[3] = 1;
            break;
        case 50:
            dsigmaydp[4] = 1;
            break;
        case 51:
            dsigmaydp[5] = 1;
            break;
        case 52:
            dsigmaydp[6] = 1;
            break;
        case 53:
            dsigmaydp[7] = 1;
            break;
        case 54:
            dsigmaydp[8] = 1;
            break;
        case 55:
            dsigmaydp[9] = 1;
            break;
        case 56:
            dsigmaydp[10] = 1;
            break;
        case 57:
            dsigmaydp[11] = 1;
            break;
        case 58:
            dsigmaydp[12] = 1;
            break;
        case 59:
            dsigmaydp[13] = 1;
            break;
        case 60:
            dsigmaydp[14] = 1;
            break;
        case 61:
            dsigmaydp[15] = 1;
            break;
        case 62:
            dsigmaydp[16] = 1;
            break;
        case 63:
            dsigmaydp[17] = 1;
            break;
        case 64:
            dsigmaydp[18] = 1;
            break;
        case 65:
            dsigmaydp[19] = 1;
            break;
        case 66:
            dsigmaydp[20] = 1;
            break;
        case 67:
            dsigmaydp[21] = 1;
            break;
        case 68:
            dsigmaydp[22] = 1;
            break;
    }
}

} // namespace model_Isensee_JCB2018
} // namespace amici
