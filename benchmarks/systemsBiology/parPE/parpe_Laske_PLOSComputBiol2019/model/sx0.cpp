#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void sx0_Laske_PLOSComputBiol2019(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 29:
            sx0[0] = 1;
            break;
    }
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
