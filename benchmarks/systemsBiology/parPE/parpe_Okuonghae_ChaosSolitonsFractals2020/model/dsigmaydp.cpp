#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void dsigmaydp_Okuonghae_ChaosSolitonsFractals2020(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 14:
            dsigmaydp[0] = 1;
            break;
        case 15:
            dsigmaydp[1] = 1;
            break;
    }
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
