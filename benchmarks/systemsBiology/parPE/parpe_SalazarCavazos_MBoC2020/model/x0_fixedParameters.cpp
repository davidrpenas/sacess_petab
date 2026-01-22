#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void x0_fixedParameters_SalazarCavazos_MBoC2020(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = EGFconc*NA*Vextra;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
        x0_fixedParameters[1] = EGFR_total;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
        x0_fixedParameters[2] = GRB2_f*GRB2_total__FREE*f;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 3) != reinitialization_state_idxs.cend())
        x0_fixedParameters[3] = SHC1_total__FREE*f;
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
