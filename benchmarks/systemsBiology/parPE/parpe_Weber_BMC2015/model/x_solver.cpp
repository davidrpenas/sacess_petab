#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Weber_BMC2015 {

void x_solver_Weber_BMC2015(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = PKD;
    x_solver[1] = PKDDAGa;
    x_solver[2] = PI4K3B;
    x_solver[3] = PI4K3Ba;
    x_solver[4] = CERTERa;
    x_solver[5] = CERT;
    x_solver[6] = CERTTGNa;
}

} // namespace model_Weber_BMC2015
} // namespace amici
