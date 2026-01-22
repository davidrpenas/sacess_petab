#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Rahman_MBS2016 {

void x_rdata_Rahman_MBS2016(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = susceptible;
    x_rdata[1] = infected_normal;
    x_rdata[2] = infected_moderate;
    x_rdata[3] = infected_weak;
    x_rdata[4] = treated_normal;
    x_rdata[5] = treated_moderate;
    x_rdata[6] = treated_weak;
}

} // namespace model_Rahman_MBS2016
} // namespace amici
