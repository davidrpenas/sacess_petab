#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Giordano_Nature2020 {

void x_solver_Giordano_Nature2020(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Susceptible;
    x_solver[1] = Infected;
    x_solver[2] = Diagnosed;
    x_solver[3] = Ailing;
    x_solver[4] = Recognized;
    x_solver[5] = Threatened;
    x_solver[6] = Healed;
    x_solver[7] = Extinct;
    x_solver[8] = CumulativeDiagnosed;
    x_solver[9] = DiagnosedHealed;
}

} // namespace model_Giordano_Nature2020
} // namespace amici
