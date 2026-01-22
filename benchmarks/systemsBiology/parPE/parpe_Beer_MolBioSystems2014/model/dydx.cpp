#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void dydx_Beer_MolBioSystems2014(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[7] = 1;
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
