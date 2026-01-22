#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

void x_solver_Armistead_CellDeathDis2024(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Cer;
    x_solver[1] = Sphingo;
    x_solver[2] = S1P;
    x_solver[3] = Sphinga;
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
