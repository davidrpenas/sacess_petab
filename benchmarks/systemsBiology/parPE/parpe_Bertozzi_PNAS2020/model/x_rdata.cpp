#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void x_rdata_Bertozzi_PNAS2020(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = I_;
    x_rdata[1] = R_;
    x_rdata[2] = S_;
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
