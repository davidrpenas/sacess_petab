#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Rahman_MBS2016 {

void dsigmaydp_Rahman_MBS2016(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 9:
            dsigmaydp[0] = 1;
            break;
    }
}

} // namespace model_Rahman_MBS2016
} // namespace amici
