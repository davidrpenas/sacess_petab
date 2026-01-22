#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Weber_BMC2015 {

void x_rdata_Weber_BMC2015(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = PKD;
    x_rdata[1] = PKDDAGa;
    x_rdata[2] = PI4K3B;
    x_rdata[3] = PI4K3Ba;
    x_rdata[4] = CERTERa;
    x_rdata[5] = CERT;
    x_rdata[6] = CERTTGNa;
}

} // namespace model_Weber_BMC2015
} // namespace amici
