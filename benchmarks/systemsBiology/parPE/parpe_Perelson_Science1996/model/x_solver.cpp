#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Perelson_Science1996 {

void x_solver_Perelson_Science1996(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = Tstar;
    x_solver[1] = V;
    x_solver[2] = Vin;
    x_solver[3] = Vni;
}

} // namespace model_Perelson_Science1996
} // namespace amici
