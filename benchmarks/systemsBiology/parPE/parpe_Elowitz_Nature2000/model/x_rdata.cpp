#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void x_rdata_Elowitz_Nature2000(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = X_protein;
    x_rdata[1] = Y_protein;
    x_rdata[2] = Z_protein;
    x_rdata[3] = GFP;
    x_rdata[4] = X_mRNA;
    x_rdata[5] = Y_mRNA;
    x_rdata[6] = Z_mRNA;
    x_rdata[7] = GFP_mRNA;
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
