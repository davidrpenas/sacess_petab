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
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void dydx_Okuonghae_ChaosSolitonsFractals2020(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[9] = 1;
    dydx[14] = 1;
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
