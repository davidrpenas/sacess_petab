#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Fujita_SciSignal2010 {

void x0_Fujita_SciSignal2010(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_EGFR;
    x0[3] = init_AKT;
    x0[5] = init_S6;
}

} // namespace model_Fujita_SciSignal2010
} // namespace amici
