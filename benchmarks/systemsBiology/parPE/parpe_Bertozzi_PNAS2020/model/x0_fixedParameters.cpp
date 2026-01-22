#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

void x0_fixedParameters_Bertozzi_PNAS2020(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
        x0_fixedParameters[2] = -I0_ + N_;
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
