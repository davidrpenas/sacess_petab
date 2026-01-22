#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void x0_Smith_BMCSystBiol2013(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[1] = initial_Ins_sim;
    x0[2] = 90000.0;
    x0[6] = 10000.0;
    x0[9] = 100000.0;
    x0[11] = 1000.0;
    x0[12] = 200000.0;
    x0[13] = 100000.0;
    x0[15] = 100000.0;
    x0[17] = 20000.0;
    x0[19] = 96000.0;
    x0[20] = 4000.0;
    x0[21] = 100000.0;
    x0[22] = 500000.0;
    x0[23] = initial_extracellular_ROS_sim;
    x0[27] = 100.0;
    x0[29] = 41700.0;
    x0[30] = 100.0;
    x0[33] = 50.0;
    x0[34] = 100000.0;
    x0[38] = 10000.0;
    x0[41] = 16000.0;
    x0[42] = 2000.0;
    x0[43] = 100000.0;
    x0[47] = 1000.0;
    x0[99] = initial_E2F1_sim;
    x0[103] = 10.0;
    x0[104] = 1000.0;
    x0[106] = 1000.0;
    x0[107] = 1000.0;
    x0[108] = 1000.0;
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
