#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void x_rdata_Blasi_CellSystems2016(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = x_0ac;
    x_rdata[1] = x_k8;
    x_rdata[2] = x_k5;
    x_rdata[3] = x_k12;
    x_rdata[4] = x_k16;
    x_rdata[5] = x_k5k8;
    x_rdata[6] = x_k5k12;
    x_rdata[7] = x_k5k16;
    x_rdata[8] = x_k8k12;
    x_rdata[9] = x_k8k16;
    x_rdata[10] = x_k12k16;
    x_rdata[11] = x_k5k8k12;
    x_rdata[12] = x_k5k8k16;
    x_rdata[13] = x_k5k12k16;
    x_rdata[14] = x_k8k12k16;
    x_rdata[15] = x_4ac;
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
