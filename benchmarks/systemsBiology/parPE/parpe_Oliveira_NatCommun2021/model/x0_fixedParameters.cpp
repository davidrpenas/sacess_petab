#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void x0_fixedParameters_Oliveira_NatCommun2021(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = exposed_init_concentration*population;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
        x0_fixedParameters[1] = population*(-asymptomatic_init_concentration - exposed_init_concentration - symptomatic_init_concentration + 1);
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
        x0_fixedParameters[2] = asymptomatic_init_concentration*population;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 3) != reinitialization_state_idxs.cend())
        x0_fixedParameters[3] = population*symptomatic_init_concentration;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 8) != reinitialization_state_idxs.cend())
        x0_fixedParameters[8] = population*symptomatic_init_concentration;
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
