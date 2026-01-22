#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void dydp_Smith_BMCSystBiol2013(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 16:
            dydp[0] = IRS1_TyrP_PI3K;
            break;
        case 17:
            dydp[1] = cellsurface_GLUT4;
            break;
        case 18:
            dydp[2] = InR_bound;
            break;
        case 19:
            dydp[3] = IRS1_TyrP;
            break;
        case 20:
            dydp[4] = PTP1B/PTP1B_plus_PTP1B_ox;
            break;
        case 21:
            dydp[5] = cellsurface_GLUT4;
            break;
        case 22:
            dydp[6] = cellsurface_GLUT4;
            break;
        case 23:
            dydp[7] = cytoplasm_SOD2;
            break;
        case 24:
            dydp[8] = Foxo1_all;
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
