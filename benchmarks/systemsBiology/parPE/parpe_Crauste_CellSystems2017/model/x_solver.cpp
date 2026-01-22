#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void x_solver_Crauste_CellSystems2017(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Naive;
    x_solver[1] = EarlyEffector;
    x_solver[2] = LateEffector;
    x_solver[3] = Memory;
    x_solver[4] = Pathogen;
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
