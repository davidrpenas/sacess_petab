#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

void x0_Armistead_CellDeathDis2024(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = initial_Cer_sim;
    x0[1] = initial_Sphingo_sim;
    x0[2] = initial_S1P_sim;
    x0[3] = initial_Sphinga_sim;
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
