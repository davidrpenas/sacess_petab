#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void y_Crauste_CellSystems2017(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = EarlyEffector;
    y[1] = EarlyEffector;
    y[2] = EarlyEffector;
    y[3] = EarlyEffector;
    y[4] = EarlyEffector;
    y[5] = EarlyEffector;
    y[6] = EarlyEffector;
    y[7] = LateEffector;
    y[8] = LateEffector;
    y[9] = LateEffector;
    y[10] = LateEffector;
    y[11] = LateEffector;
    y[12] = LateEffector;
    y[13] = LateEffector;
    y[14] = Memory;
    y[15] = Memory;
    y[16] = Memory;
    y[17] = Memory;
    y[18] = Memory;
    y[19] = Naive;
    y[20] = Naive;
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
