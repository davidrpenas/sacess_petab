#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dydx_Bruno_JExpBot2016(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[28] = 1;
    dydx[29] = 1;
    dydx[30] = 1;
    dydx[31] = 1;
    dydx[32] = 1;
    dydx[33] = 1;
    dydx[34] = 1;
    dydx[35] = 1;
    dydx[36] = 1;
    dydx[37] = 1;
    dydx[38] = 1;
    dydx[39] = 1;
    dydx[40] = 1;
    dydx[41] = 1;
    dydx[119] = 1;
    dydx[120] = 1;
    dydx[121] = 1;
    dydx[122] = 1;
    dydx[123] = 1;
    dydx[124] = 1;
    dydx[125] = 1;
    dydx[154] = 1;
    dydx[155] = 1;
    dydx[156] = 1;
    dydx[157] = 1;
    dydx[158] = 1;
    dydx[159] = 1;
    dydx[160] = 1;
    dydx[161] = 1;
    dydx[162] = 1;
    dydx[163] = 1;
    dydx[164] = 1;
    dydx[165] = 1;
    dydx[166] = 1;
    dydx[167] = 1;
    dydx[168] = 1;
    dydx[169] = 1;
    dydx[170] = 1;
    dydx[171] = 1;
    dydx[172] = 1;
    dydx[173] = 1;
    dydx[174] = 1;
    dydx[175] = 1;
    dydx[176] = 1;
    dydx[177] = 1;
    dydx[178] = 1;
    dydx[179] = 1;
    dydx[180] = 1;
    dydx[181] = 1;
    dydx[357] = 1;
    dydx[358] = 1;
    dydx[359] = 1;
    dydx[360] = 1;
    dydx[361] = 1;
    dydx[362] = 1;
    dydx[363] = 1;
    dydx[364] = 1;
    dydx[365] = 1;
    dydx[366] = 1;
    dydx[367] = 1;
    dydx[368] = 1;
    dydx[369] = 1;
    dydx[370] = 1;
    dydx[371] = 1;
    dydx[372] = 1;
    dydx[373] = 1;
    dydx[374] = 1;
    dydx[375] = 1;
    dydx[376] = 1;
    dydx[377] = 1;
    dydx[532] = 1;
    dydx[533] = 1;
    dydx[534] = 1;
    dydx[535] = 1;
    dydx[536] = 1;
    dydx[537] = 1;
    dydx[538] = 1;
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
