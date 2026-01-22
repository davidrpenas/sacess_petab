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
namespace model_Armistead_CellDeathDis2024 {

void xdot_Armistead_CellDeathDis2024(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dCerdt = flux_R2 - flux_R3 + flux_R4 - flux_R5;  // xdot[0]
    dSphingodt = flux_R3 - flux_R4 - flux_R6 + flux_R7;  // xdot[1]
    dS1Pdt = flux_R6 - flux_R7 - flux_R8;  // xdot[2]
    dSphingadt = flux_R1 - flux_R2;  // xdot[3]
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
