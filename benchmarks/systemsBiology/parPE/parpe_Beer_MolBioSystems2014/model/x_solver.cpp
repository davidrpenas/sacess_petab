#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void x_solver_Beer_MolBioSystems2014(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Bac;
    x_solver[1] = Glu;
    x_solver[2] = cGlu;
    x_solver[3] = Ind;
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
