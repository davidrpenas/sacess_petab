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
namespace model_Brannmark_JBC2010 {

void dydx_Brannmark_JBC2010(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[6] = observableParameter1_IR1_P;
    dydx[9] = observableParameter1_IR1_P;
    dydx[19] = observableParameter1_IRS1_P;
    dydx[20] = observableParameter1_IRS1_P_DosR;
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
