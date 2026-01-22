#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Chen_MSB2009 {

void dsigmaydp_Chen_MSB2009(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 235:
            dsigmaydp[0] = 1;
            break;
        case 236:
            dsigmaydp[1] = 1;
            break;
        case 237:
            dsigmaydp[2] = 1;
            break;
        case 238:
            dsigmaydp[3] = 1;
            break;
        case 239:
            dsigmaydp[4] = 1;
            break;
        case 240:
            dsigmaydp[5] = 1;
            break;
        case 241:
            dsigmaydp[6] = 1;
            break;
        case 242:
            dsigmaydp[7] = 1;
            break;
        case 243:
            dsigmaydp[8] = 1;
            break;
        case 244:
            dsigmaydp[9] = 1;
            break;
        case 245:
            dsigmaydp[10] = 1;
            break;
        case 246:
            dsigmaydp[11] = 1;
            break;
        case 247:
            dsigmaydp[12] = 1;
            break;
        case 248:
            dsigmaydp[13] = 1;
            break;
        case 249:
            dsigmaydp[14] = 1;
            break;
        case 250:
            dsigmaydp[15] = 1;
            break;
        case 251:
            dsigmaydp[16] = 1;
            break;
        case 252:
            dsigmaydp[17] = 1;
            break;
        case 253:
            dsigmaydp[18] = 1;
            break;
        case 254:
            dsigmaydp[19] = 1;
            break;
        case 255:
            dsigmaydp[20] = 1;
            break;
        case 256:
            dsigmaydp[21] = 1;
            break;
        case 257:
            dsigmaydp[22] = 1;
            break;
        case 258:
            dsigmaydp[23] = 1;
            break;
        case 259:
            dsigmaydp[24] = 1;
            break;
        case 260:
            dsigmaydp[25] = 1;
            break;
        case 261:
            dsigmaydp[26] = 1;
            break;
        case 262:
            dsigmaydp[27] = 1;
            break;
        case 263:
            dsigmaydp[28] = 1;
            break;
        case 264:
            dsigmaydp[29] = 1;
            break;
        case 265:
            dsigmaydp[30] = 1;
            break;
        case 266:
            dsigmaydp[31] = 1;
            break;
        case 267:
            dsigmaydp[32] = 1;
            break;
        case 268:
            dsigmaydp[33] = 1;
            break;
        case 269:
            dsigmaydp[34] = 1;
            break;
        case 270:
            dsigmaydp[35] = 1;
            break;
        case 271:
            dsigmaydp[36] = 1;
            break;
        case 272:
            dsigmaydp[37] = 1;
            break;
        case 273:
            dsigmaydp[38] = 1;
            break;
        case 274:
            dsigmaydp[39] = 1;
            break;
        case 275:
            dsigmaydp[40] = 1;
            break;
        case 276:
            dsigmaydp[41] = 1;
            break;
        case 277:
            dsigmaydp[42] = 1;
            break;
        case 278:
            dsigmaydp[43] = 1;
            break;
        case 279:
            dsigmaydp[44] = 1;
            break;
        case 280:
            dsigmaydp[45] = 1;
            break;
        case 281:
            dsigmaydp[46] = 1;
            break;
        case 282:
            dsigmaydp[47] = 1;
            break;
        case 283:
            dsigmaydp[48] = 1;
            break;
        case 284:
            dsigmaydp[49] = 1;
            break;
        case 285:
            dsigmaydp[50] = 1;
            break;
        case 286:
            dsigmaydp[51] = 1;
            break;
        case 287:
            dsigmaydp[52] = 1;
            break;
        case 288:
            dsigmaydp[53] = 1;
            break;
        case 289:
            dsigmaydp[54] = 1;
            break;
        case 290:
            dsigmaydp[55] = 1;
            break;
        case 291:
            dsigmaydp[56] = 1;
            break;
        case 292:
            dsigmaydp[57] = 1;
            break;
        case 293:
            dsigmaydp[58] = 1;
            break;
        case 294:
            dsigmaydp[59] = 1;
            break;
        case 295:
            dsigmaydp[60] = 1;
            break;
        case 296:
            dsigmaydp[61] = 1;
            break;
        case 297:
            dsigmaydp[62] = 1;
            break;
        case 298:
            dsigmaydp[63] = 1;
            break;
        case 299:
            dsigmaydp[64] = 1;
            break;
        case 300:
            dsigmaydp[65] = 1;
            break;
        case 301:
            dsigmaydp[66] = 1;
            break;
        case 302:
            dsigmaydp[67] = 1;
            break;
        case 303:
            dsigmaydp[68] = 1;
            break;
        case 304:
            dsigmaydp[69] = 1;
            break;
        case 305:
            dsigmaydp[70] = 1;
            break;
        case 306:
            dsigmaydp[71] = 1;
            break;
        case 307:
            dsigmaydp[72] = 1;
            break;
        case 308:
            dsigmaydp[73] = 1;
            break;
        case 309:
            dsigmaydp[74] = 1;
            break;
        case 310:
            dsigmaydp[75] = 1;
            break;
        case 311:
            dsigmaydp[76] = 1;
            break;
        case 312:
            dsigmaydp[77] = 1;
            break;
        case 313:
            dsigmaydp[78] = 1;
            break;
        case 314:
            dsigmaydp[79] = 1;
            break;
        case 315:
            dsigmaydp[80] = 1;
            break;
        case 316:
            dsigmaydp[81] = 1;
            break;
        case 317:
            dsigmaydp[82] = 1;
            break;
    }
}

} // namespace model_Chen_MSB2009
} // namespace amici
