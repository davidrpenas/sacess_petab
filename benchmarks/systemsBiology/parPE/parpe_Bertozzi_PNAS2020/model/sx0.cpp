#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void sx0_Bertozzi_PNAS2020(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 1:
            sx0[0] = 1;
            sx0[2] = -1;
            break;
    }
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
