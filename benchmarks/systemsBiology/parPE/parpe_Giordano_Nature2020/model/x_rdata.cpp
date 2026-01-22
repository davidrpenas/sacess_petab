#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void x_rdata_Giordano_Nature2020(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Susceptible;
    x_rdata[1] = Infected;
    x_rdata[2] = Diagnosed;
    x_rdata[3] = Ailing;
    x_rdata[4] = Recognized;
    x_rdata[5] = Threatened;
    x_rdata[6] = Healed;
    x_rdata[7] = Extinct;
    x_rdata[8] = CumulativeDiagnosed;
    x_rdata[9] = DiagnosedHealed;
}

} // namespace model_Giordano_Nature2020
} // namespace amici
