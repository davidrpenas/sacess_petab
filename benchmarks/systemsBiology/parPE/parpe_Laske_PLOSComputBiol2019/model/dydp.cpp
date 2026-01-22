#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void dydp_Laske_PLOSComputBiol2019(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 6:
            dydp[1] = 1;
            break;
        case 7:
            dydp[4] = 1;
            break;
        case 8:
            dydp[5] = 1;
            break;
        case 9:
            dydp[6] = 1;
            break;
        case 10:
            dydp[7] = 1;
            break;
        case 11:
            dydp[8] = 1;
            break;
        case 12:
            dydp[9] = 1;
            break;
        case 13:
            dydp[10] = 1;
            break;
        case 14:
            dydp[11] = 1;
            break;
        case 15:
            dydp[12] = 1;
            break;
    }
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
