#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void x0_Crauste_CellSystems2017(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 8090.0;
    x0[4] = 1.0;
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
