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
namespace model_Sneyd_PNAS2002 {

void dydx_Sneyd_PNAS2002(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 0.29160000000000003*std::pow(IPR_A + 0.11111111111111112*IPR_O, 3);
    dydx[4] = 2.6244000000000001*std::pow(IPR_A + 0.11111111111111112*IPR_O, 3);
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
