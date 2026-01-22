#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void dsigmaydp_Raia_CancerResearch2011(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 23:
            dsigmaydp[0] = 1;
            break;
        case 24:
            dsigmaydp[0] = observable_CD274mRNA;
            break;
        case 25:
            dsigmaydp[1] = 1;
            break;
        case 26:
            dsigmaydp[1] = observable_IL13_cell;
            break;
        case 27:
            dsigmaydp[2] = 1;
            break;
        case 28:
            dsigmaydp[2] = observable_RecSurf;
            break;
        case 29:
            dsigmaydp[3] = 1;
            break;
        case 30:
            dsigmaydp[3] = observable_SOCS3;
            break;
        case 31:
            dsigmaydp[4] = 1;
            break;
        case 32:
            dsigmaydp[4] = observable_SOCS3mRNA;
            break;
        case 33:
            dsigmaydp[5] = 1;
            break;
        case 34:
            dsigmaydp[5] = observable_pIL4Ra;
            break;
        case 35:
            dsigmaydp[6] = 1;
            break;
        case 36:
            dsigmaydp[6] = observable_pJAK2;
            break;
        case 37:
            dsigmaydp[7] = 1;
            break;
        case 38:
            dsigmaydp[7] = observable_pSTAT5;
            break;
    }
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
