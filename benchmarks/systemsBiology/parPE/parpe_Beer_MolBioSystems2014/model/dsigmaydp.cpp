#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void dsigmaydp_Beer_MolBioSystems2014(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 7:
            dsigmaydp[0] = 1;
            break;
        case 8:
            dsigmaydp[1] = 1;
            break;
    }
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
