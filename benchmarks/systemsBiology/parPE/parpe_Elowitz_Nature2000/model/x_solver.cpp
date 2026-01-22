#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void x_solver_Elowitz_Nature2000(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = X_protein;
    x_solver[1] = Y_protein;
    x_solver[2] = Z_protein;
    x_solver[3] = GFP;
    x_solver[4] = X_mRNA;
    x_solver[5] = Y_mRNA;
    x_solver[6] = Z_mRNA;
    x_solver[7] = GFP_mRNA;
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
