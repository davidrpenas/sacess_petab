#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void sx0_Smith_BMCSystBiol2013(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 25:
            sx0[23] = 1;
            break;
        case 26:
            sx0[1] = 1;
            break;
        case 27:
            sx0[99] = 1;
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
