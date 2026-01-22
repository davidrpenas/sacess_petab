#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

void sigmay_Armistead_CellDeathDis2024(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_Sphinga_obs = Sphinga_obs*noiseParameter1_Sphinga_obs;  // sigmay[0]
    sigma_Cer_obs = Cer_obs*noiseParameter1_Cer_obs;  // sigmay[1]
    sigma_Sphingo_obs = Sphingo_obs*noiseParameter1_Sphingo_obs;  // sigmay[2]
    sigma_S1P_obs = S1P_obs*noiseParameter1_S1P_obs;  // sigmay[3]
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
