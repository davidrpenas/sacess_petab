#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void sx0_Okuonghae_ChaosSolitonsFractals2020(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 11:
            sx0[1] = 1;
            break;
        case 12:
            sx0[2] = 1;
            break;
        case 13:
            sx0[3] = 1;
            break;
    }
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
