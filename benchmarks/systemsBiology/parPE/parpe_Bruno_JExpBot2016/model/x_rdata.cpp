#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void x_rdata_Bruno_JExpBot2016(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = bcar;
    x_rdata[1] = bcry;
    x_rdata[2] = b10;
    x_rdata[3] = bio;
    x_rdata[4] = ohb10;
    x_rdata[5] = ohbio;
    x_rdata[6] = zea;
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
