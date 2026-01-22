#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void w_Bertozzi_PNAS2020(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        beta_N = R0_*gamma_/N_;  // w[0]
    }

    // dynamic expressions
    flux_infection = I_*S_*beta_N;  // w[1]
    flux_recovery = I_*gamma_;  // w[2]
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
