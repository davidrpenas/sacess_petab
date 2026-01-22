#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void x_solver_Blasi_CellSystems2016(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = x_0ac;
    x_solver[1] = x_k8;
    x_solver[2] = x_k5;
    x_solver[3] = x_k12;
    x_solver[4] = x_k16;
    x_solver[5] = x_k5k8;
    x_solver[6] = x_k5k12;
    x_solver[7] = x_k5k16;
    x_solver[8] = x_k8k12;
    x_solver[9] = x_k8k16;
    x_solver[10] = x_k12k16;
    x_solver[11] = x_k5k8k12;
    x_solver[12] = x_k5k8k16;
    x_solver[13] = x_k5k12k16;
    x_solver[14] = x_k8k12k16;
    x_solver[15] = x_4ac;
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
