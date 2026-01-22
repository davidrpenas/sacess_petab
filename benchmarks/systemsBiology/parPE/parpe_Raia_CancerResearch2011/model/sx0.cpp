#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void sx0_Raia_CancerResearch2011(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 13:
            sx0[1] = 1;
            break;
    }
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
