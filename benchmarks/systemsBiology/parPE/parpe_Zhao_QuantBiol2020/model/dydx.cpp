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
namespace model_Zhao_QuantBiol2020 {

void dydx_Zhao_QuantBiol2020(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[1] = 1;
    dydx[7] = 1;
    dydx[9] = 1;
    dydx[13] = 1;
    dydx[14] = 1;
    dydx[15] = 1;
    dydx[19] = 1;
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
