#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void x_rdata_Zhao_QuantBiol2020(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Susceptible;
    x_rdata[1] = Unquarantined_Infected;
    x_rdata[2] = Quarantined_Infected;
    x_rdata[3] = Confirmed_Infected;
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
