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
#include "dwdx.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void dydx_Giordano_Nature2020(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[17] = 1;
    dydx[18] = 1;
    dydx[31] = 1;
    dydx[33] = 1;
    dydx[38] = 1;
    dydx[41] = 1;
    dydx[50] = 1;
    dydx[56] = 1;
    dydx[65] = 1;
}

} // namespace model_Giordano_Nature2020
} // namespace amici
