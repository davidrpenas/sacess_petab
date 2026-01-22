#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Sneyd_PNAS2002 {

void x_rdata_Sneyd_PNAS2002(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = IPR_O;
    x_rdata[1] = IPR_R;
    x_rdata[2] = IPR_I1;
    x_rdata[3] = IPR_S;
    x_rdata[4] = IPR_A;
    x_rdata[5] = IPR_I2;
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
