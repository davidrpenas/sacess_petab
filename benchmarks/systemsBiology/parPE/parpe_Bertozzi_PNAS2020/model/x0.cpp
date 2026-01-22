#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void x0_Bertozzi_PNAS2020(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = I0_;
    x0[2] = -I0_ + N_;
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
