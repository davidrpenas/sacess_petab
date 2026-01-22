#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void xdot_Bertozzi_PNAS2020(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dI_dt = flux_infection - flux_recovery;  // xdot[0]
    dR_dt = flux_recovery;  // xdot[1]
    dS_dt = -flux_infection;  // xdot[2]
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
