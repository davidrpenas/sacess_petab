#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void x0_fixedParameters_Zhao_QuantBiol2020(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = Total_Pop_China*Trigger_China + Total_Pop_Hubei*Trigger_Hubei + Total_Pop_Wuhan*Trigger_Wuhan;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
        x0_fixedParameters[1] = Trigger_China*(291*Trigger_Stage_I + 2800*Trigger_Stage_II) + Trigger_Hubei*(270*Trigger_Stage_I + 5700*Trigger_Stage_II) + Trigger_Wuhan*(258*Trigger_Stage_I + 15270*Trigger_Stage_II + 4000*Trigger_Stage_III);
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
        x0_fixedParameters[2] = 2000*Trigger_China*Trigger_Stage_II + 1500*Trigger_Hubei*Trigger_Stage_II + 5000*Trigger_Stage_III*Trigger_Wuhan;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 3) != reinitialization_state_idxs.cend())
        x0_fixedParameters[3] = 4000*Trigger_China*Trigger_Stage_II + 1600*Trigger_Hubei*Trigger_Stage_II + Trigger_Wuhan*(258*Trigger_Stage_I + 2000*Trigger_Stage_II + 36000*Trigger_Stage_III);
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
