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
namespace model_Fujita_SciSignal2010 {

void dydp_Fujita_SciSignal2010(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 16:
            dydp[0] = pAkt + pAkt_S6;
            break;
        case 17:
            dydp[1] = pAkt + pAkt_S6;
            break;
        case 18:
            dydp[2] = pAkt + pAkt_S6;
            break;
        case 19:
            dydp[3] = pAkt + pAkt_S6;
            break;
        case 20:
            dydp[4] = pAkt + pAkt_S6;
            break;
        case 21:
            dydp[5] = pAkt + pAkt_S6;
            break;
        case 22:
            dydp[6] = pAkt + pAkt_S6;
            break;
        case 23:
            dydp[7] = pAkt + pAkt_S6;
            break;
        case 24:
            dydp[8] = pAkt + pAkt_S6;
            break;
        case 25:
            dydp[9] = pAkt + pAkt_S6;
            break;
        case 26:
            dydp[10] = pAkt + pAkt_S6;
            break;
        case 27:
            dydp[11] = pAkt + pAkt_S6;
            break;
        case 28:
            dydp[12] = pAkt + pAkt_S6;
            break;
        case 29:
            dydp[13] = pAkt + pAkt_S6;
            break;
        case 30:
            dydp[14] = pAkt + pAkt_S6;
            break;
        case 31:
            dydp[15] = pAkt + pAkt_S6;
            break;
        case 32:
            dydp[16] = pAkt + pAkt_S6;
            break;
        case 33:
            dydp[17] = pAkt + pAkt_S6;
            break;
        case 34:
            dydp[18] = pAkt + pAkt_S6;
            break;
        case 35:
            dydp[19] = pAkt + pAkt_S6;
            break;
        case 36:
            dydp[20] = pEGFR + pEGFR_Akt;
            break;
        case 37:
            dydp[21] = pEGFR + pEGFR_Akt;
            break;
        case 38:
            dydp[22] = pEGFR + pEGFR_Akt;
            break;
        case 39:
            dydp[23] = pEGFR + pEGFR_Akt;
            break;
        case 40:
            dydp[24] = pEGFR + pEGFR_Akt;
            break;
        case 41:
            dydp[25] = pEGFR + pEGFR_Akt;
            break;
        case 42:
            dydp[26] = pEGFR + pEGFR_Akt;
            break;
        case 43:
            dydp[27] = pEGFR + pEGFR_Akt;
            break;
        case 44:
            dydp[28] = pEGFR + pEGFR_Akt;
            break;
        case 45:
            dydp[29] = pEGFR + pEGFR_Akt;
            break;
        case 46:
            dydp[30] = pEGFR + pEGFR_Akt;
            break;
        case 47:
            dydp[31] = pEGFR + pEGFR_Akt;
            break;
        case 48:
            dydp[32] = pEGFR + pEGFR_Akt;
            break;
        case 49:
            dydp[33] = pEGFR + pEGFR_Akt;
            break;
        case 50:
            dydp[34] = pEGFR + pEGFR_Akt;
            break;
        case 51:
            dydp[35] = pEGFR + pEGFR_Akt;
            break;
        case 52:
            dydp[36] = pEGFR + pEGFR_Akt;
            break;
        case 53:
            dydp[37] = pEGFR + pEGFR_Akt;
            break;
        case 54:
            dydp[38] = pEGFR + pEGFR_Akt;
            break;
        case 55:
            dydp[39] = pEGFR + pEGFR_Akt;
            break;
        case 56:
            dydp[40] = pEGFR + pEGFR_Akt;
            break;
        case 57:
            dydp[41] = pEGFR + pEGFR_Akt;
            break;
        case 58:
            dydp[42] = pEGFR + pEGFR_Akt;
            break;
        case 59:
            dydp[43] = pS6;
            break;
        case 60:
            dydp[44] = pS6;
            break;
        case 61:
            dydp[45] = pS6;
            break;
        case 62:
            dydp[46] = pS6;
            break;
        case 63:
            dydp[47] = pS6;
            break;
        case 64:
            dydp[48] = pS6;
            break;
        case 65:
            dydp[49] = pS6;
            break;
        case 66:
            dydp[50] = pS6;
            break;
        case 67:
            dydp[51] = pS6;
            break;
        case 68:
            dydp[52] = pS6;
            break;
        case 69:
            dydp[53] = pS6;
            break;
        case 70:
            dydp[54] = pS6;
            break;
        case 71:
            dydp[55] = pS6;
            break;
        case 72:
            dydp[56] = pS6;
            break;
        case 73:
            dydp[57] = pS6;
            break;
        case 74:
            dydp[58] = pS6;
            break;
        case 75:
            dydp[59] = pS6;
            break;
        case 76:
            dydp[60] = pS6;
            break;
        case 77:
            dydp[61] = pS6;
            break;
        case 78:
            dydp[62] = pS6;
            break;
        case 79:
            dydp[63] = pS6;
            break;
        case 80:
            dydp[64] = pS6;
            break;
        case 81:
            dydp[65] = pS6;
            break;
        case 82:
            dydp[66] = pS6;
            break;
    }
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
