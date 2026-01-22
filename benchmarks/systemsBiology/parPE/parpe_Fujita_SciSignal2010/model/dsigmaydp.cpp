#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Fujita_SciSignal2010 {

void dsigmaydp_Fujita_SciSignal2010(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 83:
            dsigmaydp[0] = 1;
            break;
        case 84:
            dsigmaydp[1] = 1;
            break;
        case 85:
            dsigmaydp[2] = 1;
            break;
        case 86:
            dsigmaydp[3] = 1;
            break;
        case 87:
            dsigmaydp[4] = 1;
            break;
        case 88:
            dsigmaydp[5] = 1;
            break;
        case 89:
            dsigmaydp[6] = 1;
            break;
        case 90:
            dsigmaydp[7] = 1;
            break;
        case 91:
            dsigmaydp[8] = 1;
            break;
        case 92:
            dsigmaydp[9] = 1;
            break;
        case 93:
            dsigmaydp[10] = 1;
            break;
        case 94:
            dsigmaydp[11] = 1;
            break;
        case 95:
            dsigmaydp[12] = 1;
            break;
        case 96:
            dsigmaydp[13] = 1;
            break;
        case 97:
            dsigmaydp[14] = 1;
            break;
        case 98:
            dsigmaydp[15] = 1;
            break;
        case 99:
            dsigmaydp[16] = 1;
            break;
        case 100:
            dsigmaydp[17] = 1;
            break;
        case 101:
            dsigmaydp[18] = 1;
            break;
        case 102:
            dsigmaydp[19] = 1;
            break;
        case 103:
            dsigmaydp[20] = 1;
            break;
        case 104:
            dsigmaydp[21] = 1;
            break;
        case 105:
            dsigmaydp[22] = 1;
            break;
        case 106:
            dsigmaydp[23] = 1;
            break;
        case 107:
            dsigmaydp[24] = 1;
            break;
        case 108:
            dsigmaydp[25] = 1;
            break;
        case 109:
            dsigmaydp[26] = 1;
            break;
        case 110:
            dsigmaydp[27] = 1;
            break;
        case 111:
            dsigmaydp[28] = 1;
            break;
        case 112:
            dsigmaydp[29] = 1;
            break;
        case 113:
            dsigmaydp[30] = 1;
            break;
        case 114:
            dsigmaydp[31] = 1;
            break;
        case 115:
            dsigmaydp[32] = 1;
            break;
        case 116:
            dsigmaydp[33] = 1;
            break;
        case 117:
            dsigmaydp[34] = 1;
            break;
        case 118:
            dsigmaydp[35] = 1;
            break;
        case 119:
            dsigmaydp[36] = 1;
            break;
        case 120:
            dsigmaydp[37] = 1;
            break;
        case 121:
            dsigmaydp[38] = 1;
            break;
        case 122:
            dsigmaydp[39] = 1;
            break;
        case 123:
            dsigmaydp[40] = 1;
            break;
        case 124:
            dsigmaydp[41] = 1;
            break;
        case 125:
            dsigmaydp[42] = 1;
            break;
        case 126:
            dsigmaydp[43] = 1;
            break;
        case 127:
            dsigmaydp[44] = 1;
            break;
        case 128:
            dsigmaydp[45] = 1;
            break;
        case 129:
            dsigmaydp[46] = 1;
            break;
        case 130:
            dsigmaydp[47] = 1;
            break;
        case 131:
            dsigmaydp[48] = 1;
            break;
        case 132:
            dsigmaydp[49] = 1;
            break;
        case 133:
            dsigmaydp[50] = 1;
            break;
        case 134:
            dsigmaydp[51] = 1;
            break;
        case 135:
            dsigmaydp[52] = 1;
            break;
        case 136:
            dsigmaydp[53] = 1;
            break;
        case 137:
            dsigmaydp[54] = 1;
            break;
        case 138:
            dsigmaydp[55] = 1;
            break;
        case 139:
            dsigmaydp[56] = 1;
            break;
        case 140:
            dsigmaydp[57] = 1;
            break;
        case 141:
            dsigmaydp[58] = 1;
            break;
        case 142:
            dsigmaydp[59] = 1;
            break;
        case 143:
            dsigmaydp[60] = 1;
            break;
        case 144:
            dsigmaydp[61] = 1;
            break;
        case 145:
            dsigmaydp[62] = 1;
            break;
        case 146:
            dsigmaydp[63] = 1;
            break;
        case 147:
            dsigmaydp[64] = 1;
            break;
        case 148:
            dsigmaydp[65] = 1;
            break;
        case 149:
            dsigmaydp[66] = 1;
            break;
    }
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
