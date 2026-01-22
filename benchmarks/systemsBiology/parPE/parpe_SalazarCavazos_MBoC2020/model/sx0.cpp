#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void sx0_SalazarCavazos_MBoC2020(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 0:
            sx0[2] = GRB2_f*f;
            break;
        case 1:
            sx0[3] = f;
            break;
    }
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
