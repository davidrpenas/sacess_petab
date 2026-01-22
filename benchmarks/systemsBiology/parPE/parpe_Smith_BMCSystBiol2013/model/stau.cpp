#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "sx.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void stau_Smith_BMCSystBiol2013(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *sx, const int ip, const int ie){
    switch(ie) {
        case 1:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                case 19:
                case 20:
                case 21:
                case 22:
                case 23:
                case 24:
                case 25:
                case 26:
                case 27:
                    stau[0] = -sx11/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    break;
            }
            break;
        case 2:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                case 19:
                case 20:
                case 21:
                case 22:
                case 23:
                case 24:
                case 25:
                case 26:
                case 27:
                    stau[0] = sx11/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    break;
            }
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
