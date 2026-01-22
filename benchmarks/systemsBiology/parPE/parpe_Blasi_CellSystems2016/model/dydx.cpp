#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void dydx_Blasi_CellSystems2016(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[28] = 1;
    dydx[37] = 1;
    dydx[50] = 1;
    dydx[68] = 1;
    dydx[89] = 1;
    dydx[102] = 1;
    dydx[120] = 1;
    dydx[141] = 1;
    dydx[159] = 1;
    dydx[163] = 1;
    dydx[186] = 1;
    dydx[203] = 1;
    dydx[215] = 1;
    dydx[238] = 1;
    dydx[241] = 1;
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
