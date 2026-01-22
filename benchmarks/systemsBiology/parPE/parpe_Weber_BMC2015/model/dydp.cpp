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
namespace model_Weber_BMC2015 {

void dydp_Weber_BMC2015(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 26:
            dydp[0] = CERT/(CERT + CERTERa + CERTTGNa);
            break;
        case 27:
            dydp[2] = PI4K3Ba/(PI4K3B + PI4K3Ba);
            break;
        case 28:
            dydp[4] = PKDDAGa;
            break;
        case 29:
            dydp[5] = PKDDAGa;
            break;
        case 30:
            dydp[6] = PKDDAGa;
            break;
    }
}

} // namespace model_Weber_BMC2015
} // namespace amici
