#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void sx0_Oliveira_NatCommun2021(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 9:
            sx0[0] = population;
            sx0[1] = -population;
            break;
        case 10:
            sx0[1] = -population;
            sx0[2] = population;
            break;
        case 11:
            sx0[1] = -population;
            sx0[3] = population;
            sx0[8] = population;
            break;
    }
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
