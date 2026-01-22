#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void dsigmaydp_Elowitz_Nature2000(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 20:
            dsigmaydp[0] = 1;
            break;
    }
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
