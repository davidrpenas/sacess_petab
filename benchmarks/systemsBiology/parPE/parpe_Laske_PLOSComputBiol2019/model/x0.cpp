#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void x0_Laske_PLOSComputBiol2019(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = initial_V_ex_sim;
    x0[32] = 1000.0;
    x0[33] = 150.0;
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
