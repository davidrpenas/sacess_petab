#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void x_rdata_Okuonghae_ChaosSolitonsFractals2020(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = susceptible;
    x_rdata[1] = exposed;
    x_rdata[2] = asymptomatic;
    x_rdata[3] = symptomatic;
    x_rdata[4] = detected;
    x_rdata[5] = recovered;
    x_rdata[6] = deceased;
    x_rdata[7] = detected_cumulative;
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
