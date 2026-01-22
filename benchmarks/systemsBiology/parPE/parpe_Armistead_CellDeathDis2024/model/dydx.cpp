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
namespace model_Armistead_CellDeathDis2024 {

void dydx_Armistead_CellDeathDis2024(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[1] = 1;
    dydx[6] = 1;
    dydx[11] = 1;
    dydx[12] = 1;
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
