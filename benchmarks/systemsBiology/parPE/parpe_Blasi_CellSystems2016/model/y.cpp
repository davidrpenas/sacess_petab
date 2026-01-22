#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void y_Blasi_CellSystems2016(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = x_0ac;
    y[1] = x_4ac;
    y[2] = x_k12;
    y[3] = x_k12k16;
    y[4] = x_k16;
    y[5] = x_k5;
    y[6] = x_k5k12;
    y[7] = x_k5k12k16;
    y[8] = x_k5k16;
    y[9] = x_k5k8;
    y[10] = x_k5k8k12;
    y[11] = x_k5k8k16;
    y[12] = x_k8;
    y[13] = x_k8k12;
    y[14] = x_k8k12k16;
    y[15] = x_k8k16;
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
