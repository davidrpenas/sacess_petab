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
namespace model_Raia_CancerResearch2011 {

void dydx_Raia_CancerResearch2011(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[2] = 1;
    dydx[17] = observableParameter1_observable_IL13_cell;
    dydx[18] = 1;
    dydx[25] = observableParameter1_observable_IL13_cell;
    dydx[26] = 1;
    dydx[29] = observableParameter1_observable_pIL4Ra;
    dydx[33] = observableParameter1_observable_IL13_cell;
    dydx[37] = observableParameter1_observable_pIL4Ra;
    dydx[54] = observableParameter1_observable_pJAK2;
    dydx[71] = 1;
    dydx[76] = observableParameter1_observable_SOCS3mRNA;
    dydx[89] = observableParameter1_observable_IL13_cell;
    dydx[99] = 1;
    dydx[104] = observableParameter1_observable_CD274mRNA;
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
