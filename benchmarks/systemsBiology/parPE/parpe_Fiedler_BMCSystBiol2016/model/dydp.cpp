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
namespace model_Fiedler_BMCSystBiol2016 {

void dydp_Fiedler_BMCSystBiol2016(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 12:
            dydp[0] = pERK;
            break;
        case 13:
            dydp[1] = pERK;
            break;
        case 14:
            dydp[2] = pERK;
            break;
        case 15:
            dydp[3] = pERK;
            break;
        case 16:
            dydp[4] = pERK;
            break;
        case 17:
            dydp[5] = pERK;
            break;
        case 18:
            dydp[6] = pMEK;
            break;
        case 19:
            dydp[7] = pMEK;
            break;
        case 20:
            dydp[8] = pMEK;
            break;
        case 21:
            dydp[9] = pMEK;
            break;
        case 22:
            dydp[10] = pMEK;
            break;
        case 23:
            dydp[11] = pMEK;
            break;
    }
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
