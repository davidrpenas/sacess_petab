#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Fujita_SciSignal2010 {

void x_rdata_Fujita_SciSignal2010(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = EGFR;
    x_rdata[1] = pEGFR;
    x_rdata[2] = pEGFR_Akt;
    x_rdata[3] = Akt;
    x_rdata[4] = pAkt;
    x_rdata[5] = S6;
    x_rdata[6] = pAkt_S6;
    x_rdata[7] = pS6;
    x_rdata[8] = EGF_EGFR;
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
