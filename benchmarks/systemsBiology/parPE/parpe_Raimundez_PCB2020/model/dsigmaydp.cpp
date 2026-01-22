#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void dsigmaydp_Raimundez_PCB2020(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 136:
            dsigmaydp[0] = 1;
            break;
        case 137:
            dsigmaydp[1] = 1;
            break;
        case 138:
            dsigmaydp[2] = 1;
            break;
        case 139:
            dsigmaydp[3] = 1;
            break;
        case 140:
            dsigmaydp[4] = 1;
            break;
        case 141:
            dsigmaydp[5] = 1;
            break;
        case 142:
            dsigmaydp[6] = 1;
            break;
        case 143:
            dsigmaydp[7] = 1;
            break;
        case 144:
            dsigmaydp[8] = 1;
            break;
        case 145:
            dsigmaydp[9] = 1;
            break;
        case 146:
            dsigmaydp[10] = 1;
            break;
        case 147:
            dsigmaydp[11] = 1;
            break;
        case 148:
            dsigmaydp[12] = 1;
            break;
        case 149:
            dsigmaydp[13] = 1;
            break;
        case 150:
            dsigmaydp[14] = 1;
            break;
        case 151:
            dsigmaydp[15] = 1;
            break;
        case 152:
            dsigmaydp[16] = 1;
            break;
        case 153:
            dsigmaydp[17] = 1;
            break;
        case 154:
            dsigmaydp[18] = 1;
            break;
        case 155:
            dsigmaydp[19] = 1;
            break;
        case 156:
            dsigmaydp[20] = 1;
            break;
        case 157:
            dsigmaydp[21] = 1;
            break;
        case 158:
            dsigmaydp[22] = 1;
            break;
        case 159:
            dsigmaydp[23] = 1;
            break;
        case 160:
            dsigmaydp[24] = 1;
            break;
        case 161:
            dsigmaydp[25] = 1;
            break;
        case 162:
            dsigmaydp[26] = 1;
            break;
        case 163:
            dsigmaydp[27] = 1;
            break;
        case 164:
            dsigmaydp[28] = 1;
            break;
        case 165:
            dsigmaydp[29] = 1;
            break;
        case 166:
            dsigmaydp[30] = 1;
            break;
        case 167:
            dsigmaydp[31] = 1;
            break;
        case 168:
            dsigmaydp[32] = 1;
            break;
        case 169:
            dsigmaydp[33] = 1;
            break;
        case 170:
            dsigmaydp[34] = 1;
            break;
        case 171:
            dsigmaydp[35] = 1;
            break;
        case 172:
            dsigmaydp[36] = 1;
            break;
        case 173:
            dsigmaydp[37] = 1;
            break;
        case 174:
            dsigmaydp[38] = 1;
            break;
        case 175:
            dsigmaydp[39] = 1;
            break;
        case 176:
            dsigmaydp[40] = 1;
            break;
        case 177:
            dsigmaydp[41] = 1;
            break;
        case 178:
            dsigmaydp[42] = 1;
            break;
        case 179:
            dsigmaydp[43] = 1;
            break;
        case 180:
            dsigmaydp[44] = 1;
            break;
        case 181:
            dsigmaydp[45] = 1;
            break;
        case 182:
            dsigmaydp[46] = 1;
            break;
        case 183:
            dsigmaydp[47] = 1;
            break;
        case 184:
            dsigmaydp[48] = 1;
            break;
        case 185:
            dsigmaydp[49] = 1;
            break;
        case 186:
            dsigmaydp[50] = 1;
            break;
        case 187:
            dsigmaydp[51] = 1;
            break;
        case 188:
            dsigmaydp[52] = 1;
            break;
        case 189:
            dsigmaydp[53] = 1;
            break;
        case 190:
            dsigmaydp[54] = 1;
            break;
        case 191:
            dsigmaydp[55] = 1;
            break;
        case 192:
            dsigmaydp[56] = 1;
            break;
        case 193:
            dsigmaydp[57] = 1;
            break;
        case 194:
            dsigmaydp[58] = 1;
            break;
        case 195:
            dsigmaydp[59] = 1;
            break;
        case 196:
            dsigmaydp[60] = 1;
            break;
        case 197:
            dsigmaydp[61] = 1;
            break;
        case 198:
            dsigmaydp[62] = 1;
            break;
        case 199:
            dsigmaydp[63] = 1;
            break;
        case 200:
            dsigmaydp[64] = 1;
            break;
        case 201:
            dsigmaydp[65] = 1;
            break;
        case 202:
            dsigmaydp[66] = 1;
            break;
        case 203:
            dsigmaydp[67] = 1;
            break;
        case 204:
            dsigmaydp[68] = 1;
            break;
        case 205:
            dsigmaydp[69] = 1;
            break;
        case 206:
            dsigmaydp[70] = 1;
            break;
        case 207:
            dsigmaydp[71] = 1;
            break;
        case 208:
            dsigmaydp[72] = 1;
            break;
        case 209:
            dsigmaydp[73] = 1;
            break;
        case 210:
            dsigmaydp[74] = 1;
            break;
        case 211:
            dsigmaydp[75] = 1;
            break;
        case 212:
            dsigmaydp[76] = 1;
            break;
        case 213:
            dsigmaydp[77] = 1;
            break;
        case 214:
            dsigmaydp[78] = 1;
            break;
        case 215:
            dsigmaydp[79] = 1;
            break;
        case 216:
            dsigmaydp[80] = 1;
            break;
        case 217:
            dsigmaydp[81] = 1;
            break;
        case 218:
            dsigmaydp[82] = 1;
            break;
        case 219:
            dsigmaydp[83] = 1;
            break;
        case 220:
            dsigmaydp[84] = 1;
            break;
        case 221:
            dsigmaydp[85] = 1;
            break;
        case 222:
            dsigmaydp[86] = 1;
            break;
        case 223:
            dsigmaydp[87] = 1;
            break;
        case 224:
            dsigmaydp[88] = 1;
            break;
        case 225:
            dsigmaydp[89] = 1;
            break;
        case 226:
            dsigmaydp[90] = 1;
            break;
        case 227:
            dsigmaydp[91] = 1;
            break;
        case 228:
            dsigmaydp[92] = 1;
            break;
        case 229:
            dsigmaydp[93] = 1;
            break;
        case 230:
            dsigmaydp[94] = 1;
            break;
        case 231:
            dsigmaydp[95] = 1;
            break;
        case 232:
            dsigmaydp[96] = 1;
            break;
        case 233:
            dsigmaydp[97] = 1;
            break;
        case 234:
            dsigmaydp[98] = 1;
            break;
        case 235:
            dsigmaydp[99] = 1;
            break;
        case 236:
            dsigmaydp[100] = 1;
            break;
        case 237:
            dsigmaydp[101] = 1;
            break;
        case 238:
            dsigmaydp[102] = 1;
            break;
        case 239:
            dsigmaydp[103] = 1;
            break;
        case 240:
            dsigmaydp[104] = 1;
            break;
        case 241:
            dsigmaydp[105] = 1;
            break;
        case 242:
            dsigmaydp[106] = 1;
            break;
        case 243:
            dsigmaydp[107] = 1;
            break;
        case 244:
            dsigmaydp[108] = 1;
            break;
        case 245:
            dsigmaydp[109] = 1;
            break;
        case 246:
            dsigmaydp[110] = 1;
            break;
    }
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
