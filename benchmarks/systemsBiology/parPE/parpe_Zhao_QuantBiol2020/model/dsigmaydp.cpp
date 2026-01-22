#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void dsigmaydp_Zhao_QuantBiol2020(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 21:
            dsigmaydp[0] = 1;
            dsigmaydp[2] = 1;
            break;
        case 22:
            dsigmaydp[1] = 1;
            break;
        case 23:
            dsigmaydp[3] = 1;
            break;
        case 24:
            dsigmaydp[4] = 1;
            break;
    }
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
