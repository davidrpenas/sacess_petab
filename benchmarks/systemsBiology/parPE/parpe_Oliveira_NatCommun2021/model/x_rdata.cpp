#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void x_rdata_Oliveira_NatCommun2021(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Exposed;
    x_rdata[1] = Susceptible;
    x_rdata[2] = Asymptomatic;
    x_rdata[3] = Symptomatic;
    x_rdata[4] = Recovered;
    x_rdata[5] = Hospital;
    x_rdata[6] = ICU;
    x_rdata[7] = Deaths;
    x_rdata[8] = Cumulative_cases;
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
