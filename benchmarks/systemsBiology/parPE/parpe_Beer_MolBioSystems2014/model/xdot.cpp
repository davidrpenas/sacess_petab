#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "w.h"
#include "xdot.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void xdot_Beer_MolBioSystems2014(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dBacdt = flux_reaction1;  // xdot[0]
    dGludt = flux_reaction2;  // xdot[1]
    dcGludt = flux_reaction3;  // xdot[2]
    dInddt = flux_reaction4;  // xdot[3]
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
