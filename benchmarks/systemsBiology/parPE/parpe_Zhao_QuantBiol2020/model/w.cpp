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
namespace model_Zhao_QuantBiol2020 {

void w_Zhao_QuantBiol2020(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        R = Trigger_China*(R_Stage_II_China*Trigger_Stage_II + R_Stage_I_China*Trigger_Stage_I) + Trigger_Hubei*(R_Stage_II_Hubei*Trigger_Stage_II + R_Stage_I_Hubei*Trigger_Stage_I) + Trigger_Wuhan*(R_Stage_III_Wuhan*Trigger_Stage_III + R_Stage_II_Wuhan*Trigger_Stage_II + R_Stage_I_Wuhan*Trigger_Stage_I);  // w[1]
        Total_Pop = Total_Pop_China*Trigger_China + Total_Pop_Hubei*Trigger_Hubei + Total_Pop_Wuhan*Trigger_Wuhan;  // w[2]
        gamma_1 = Trigger_China*(Trigger_Stage_I*gamma_1_Stage_I_China + Trigger_Stage_II*gamma_1_Stage_II_China) + Trigger_Hubei*(Trigger_Stage_I*gamma_1_Stage_I_Hubei + Trigger_Stage_II*gamma_1_Stage_II_Hubei) + Trigger_Wuhan*(Trigger_Stage_I*gamma_1_Stage_I_Wuhan + Trigger_Stage_II*gamma_1_Stage_II_Wuhan + Trigger_Stage_III*gamma_1_Stage_III_Wuhan);  // w[3]
        gamma_2 = Trigger_China*(Trigger_Stage_I*gamma_2_Stage_I_China + Trigger_Stage_II*gamma_2_Stage_II_China) + Trigger_Hubei*(Trigger_Stage_I*gamma_2_Stage_I_Hubei + Trigger_Stage_II*gamma_2_Stage_II_Hubei) + Trigger_Wuhan*(Trigger_Stage_I*gamma_2_Stage_I_Wuhan + Trigger_Stage_II*gamma_2_Stage_II_Wuhan + Trigger_Stage_III*gamma_2_Stage_III_Wuhan);  // w[4]
        alpha = R*gamma_1;  // w[5]
        beta = gamma_2 + sigma*(1 - gamma_2);  // w[6]
    }

    // dynamic expressions
    Cumulative_Infected = Confirmed_Infected + Quarantined_Infected + Unquarantined_Infected;  // w[0]
    flux_Susceptible_to_Unquarantined = 1.0*Susceptible*Unquarantined_Infected*alpha/Total_Pop;  // w[7]
    flux_Unquarantined_to_Quarantined = 1.0*Unquarantined_Infected*gamma_1;  // w[8]
    flux_Quarantined_to_Confirmed = 1.0*Quarantined_Infected*beta;  // w[9]
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
