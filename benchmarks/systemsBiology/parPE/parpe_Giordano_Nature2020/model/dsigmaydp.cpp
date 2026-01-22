#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void dsigmaydp_Giordano_Nature2020(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 43:
            dsigmaydp[0] = 1;
            break;
        case 44:
            dsigmaydp[1] = 1;
            break;
        case 45:
            dsigmaydp[2] = 1;
            break;
        case 46:
            dsigmaydp[3] = 1;
            break;
        case 47:
            dsigmaydp[4] = 1;
            break;
        case 48:
            dsigmaydp[5] = 1;
            break;
        case 49:
            dsigmaydp[6] = 1;
            break;
    }
}

} // namespace model_Giordano_Nature2020
} // namespace amici
