#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Fujita_SciSignal2010 {

void x_solver_Fujita_SciSignal2010(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EGFR;
    x_solver[1] = pEGFR;
    x_solver[2] = pEGFR_Akt;
    x_solver[3] = Akt;
    x_solver[4] = pAkt;
    x_solver[5] = S6;
    x_solver[6] = pAkt_S6;
    x_solver[7] = pS6;
    x_solver[8] = EGF_EGFR;
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
