#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void x_rdata_Fiedler_BMCSystBiol2016(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = RAF;
    x_rdata[1] = pRAF;
    x_rdata[2] = MEK;
    x_rdata[3] = pMEK;
    x_rdata[4] = ERK;
    x_rdata[5] = pERK;
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
