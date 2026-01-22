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
namespace model_Crauste_CellSystems2017 {

void dydx_Crauste_CellSystems2017(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[19] = 1;
    dydx[20] = 1;
    dydx[21] = 1;
    dydx[22] = 1;
    dydx[23] = 1;
    dydx[24] = 1;
    dydx[25] = 1;
    dydx[26] = 1;
    dydx[27] = 1;
    dydx[49] = 1;
    dydx[50] = 1;
    dydx[51] = 1;
    dydx[52] = 1;
    dydx[53] = 1;
    dydx[54] = 1;
    dydx[55] = 1;
    dydx[77] = 1;
    dydx[78] = 1;
    dydx[79] = 1;
    dydx[80] = 1;
    dydx[81] = 1;
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
