#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void x0_Lucarelli_CellSystems2018(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_TGFb;
    x0[1] = init_Rec;
    x0[3] = init_S2;
    x0[4] = init_S3;
    x0[5] = init_S4;
    x0[21] = 1.0;
    x0[22] = 1.0;
    x0[23] = 1.0;
    x0[24] = 1.0;
    x0[25] = 1.0;
    x0[26] = 1.0;
    x0[27] = 1.0;
    x0[28] = 1.0;
    x0[29] = 1.0;
    x0[30] = 1.0;
    x0[31] = 1.0;
    x0[32] = 1.0;
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
