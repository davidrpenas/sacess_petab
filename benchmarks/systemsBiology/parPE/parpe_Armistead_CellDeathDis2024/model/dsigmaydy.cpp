#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

void dsigmaydy_Armistead_CellDeathDis2024(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    dsigmaydy[0] = noiseParameter1_Sphinga_obs;
    dsigmaydy[5] = noiseParameter1_Cer_obs;
    dsigmaydy[10] = noiseParameter1_Sphingo_obs;
    dsigmaydy[15] = noiseParameter1_S1P_obs;
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
