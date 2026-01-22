#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dsigmaydp_Bruno_JExpBot2016(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 17:
            dsigmaydp[0] = 1;
            break;
        case 18:
            dsigmaydp[1] = 1;
            break;
        case 19:
            dsigmaydp[2] = 1;
            break;
        case 20:
            dsigmaydp[3] = 1;
            break;
        case 21:
            dsigmaydp[4] = 1;
            break;
        case 22:
            dsigmaydp[5] = 1;
            break;
        case 23:
            dsigmaydp[6] = 1;
            break;
        case 24:
            dsigmaydp[7] = 1;
            break;
        case 25:
            dsigmaydp[8] = 1;
            break;
        case 26:
            dsigmaydp[9] = 1;
            break;
        case 27:
            dsigmaydp[10] = 1;
            break;
        case 28:
            dsigmaydp[11] = 1;
            break;
        case 29:
            dsigmaydp[12] = 1;
            break;
        case 30:
            dsigmaydp[13] = 1;
            break;
        case 31:
            dsigmaydp[14] = 1;
            break;
        case 32:
            dsigmaydp[15] = 1;
            break;
        case 33:
            dsigmaydp[16] = 1;
            break;
        case 34:
            dsigmaydp[17] = 1;
            break;
        case 35:
            dsigmaydp[18] = 1;
            break;
        case 36:
            dsigmaydp[19] = 1;
            break;
        case 37:
            dsigmaydp[20] = 1;
            break;
        case 38:
            dsigmaydp[21] = 1;
            break;
        case 39:
            dsigmaydp[22] = 1;
            break;
        case 40:
            dsigmaydp[23] = 1;
            break;
        case 41:
            dsigmaydp[24] = 1;
            break;
        case 42:
            dsigmaydp[25] = 1;
            break;
        case 43:
            dsigmaydp[26] = 1;
            break;
        case 44:
            dsigmaydp[27] = 1;
            break;
        case 45:
            dsigmaydp[28] = 1;
            break;
        case 46:
            dsigmaydp[29] = 1;
            break;
        case 47:
            dsigmaydp[30] = 1;
            break;
        case 48:
            dsigmaydp[31] = 1;
            break;
        case 49:
            dsigmaydp[32] = 1;
            break;
        case 50:
            dsigmaydp[33] = 1;
            break;
        case 51:
            dsigmaydp[34] = 1;
            break;
        case 52:
            dsigmaydp[35] = 1;
            break;
        case 53:
            dsigmaydp[36] = 1;
            break;
        case 54:
            dsigmaydp[37] = 1;
            break;
        case 55:
            dsigmaydp[38] = 1;
            break;
        case 56:
            dsigmaydp[39] = 1;
            break;
        case 57:
            dsigmaydp[40] = 1;
            break;
        case 58:
            dsigmaydp[41] = 1;
            break;
        case 59:
            dsigmaydp[42] = 1;
            break;
        case 60:
            dsigmaydp[43] = 1;
            break;
        case 61:
            dsigmaydp[44] = 1;
            break;
        case 62:
            dsigmaydp[45] = 1;
            break;
        case 63:
            dsigmaydp[46] = 1;
            break;
        case 64:
            dsigmaydp[47] = 1;
            break;
        case 65:
            dsigmaydp[48] = 1;
            break;
        case 66:
            dsigmaydp[49] = 1;
            break;
        case 67:
            dsigmaydp[50] = 1;
            break;
        case 68:
            dsigmaydp[51] = 1;
            break;
        case 69:
            dsigmaydp[52] = 1;
            break;
        case 70:
            dsigmaydp[53] = 1;
            break;
        case 71:
            dsigmaydp[54] = 1;
            break;
        case 72:
            dsigmaydp[55] = 1;
            break;
        case 73:
            dsigmaydp[56] = 1;
            break;
        case 74:
            dsigmaydp[57] = 1;
            break;
        case 75:
            dsigmaydp[58] = 1;
            break;
        case 76:
            dsigmaydp[59] = 1;
            break;
        case 77:
            dsigmaydp[60] = 1;
            break;
        case 78:
            dsigmaydp[61] = 1;
            break;
        case 79:
            dsigmaydp[62] = 1;
            break;
        case 80:
            dsigmaydp[63] = 1;
            break;
        case 81:
            dsigmaydp[64] = 1;
            break;
        case 82:
            dsigmaydp[65] = 1;
            break;
        case 83:
            dsigmaydp[66] = 1;
            break;
        case 84:
            dsigmaydp[67] = 1;
            break;
        case 85:
            dsigmaydp[68] = 1;
            break;
        case 86:
            dsigmaydp[69] = 1;
            break;
        case 87:
            dsigmaydp[70] = 1;
            break;
        case 88:
            dsigmaydp[71] = 1;
            break;
        case 89:
            dsigmaydp[72] = 1;
            break;
        case 90:
            dsigmaydp[73] = 1;
            break;
        case 91:
            dsigmaydp[74] = 1;
            break;
        case 92:
            dsigmaydp[75] = 1;
            break;
        case 93:
            dsigmaydp[76] = 1;
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
