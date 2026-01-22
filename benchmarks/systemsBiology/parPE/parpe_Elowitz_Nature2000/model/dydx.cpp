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
namespace model_Elowitz_Nature2000 {

void dydx_Elowitz_Nature2000(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[3] = observableParameter2_fluorescence;
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
