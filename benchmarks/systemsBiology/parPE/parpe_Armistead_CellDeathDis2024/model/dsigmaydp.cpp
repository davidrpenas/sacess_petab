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

void dsigmaydp_Armistead_CellDeathDis2024(realtype *dsigmaydp, const realtype t, const realtype *p, const realtype *k, const realtype *y, const int ip){
    switch(ip) {
        case 10:
            dsigmaydp[0] = Sphinga_obs;
            break;
        case 11:
            dsigmaydp[1] = Cer_obs;
            break;
        case 12:
            dsigmaydp[2] = Sphingo_obs;
            break;
        case 13:
            dsigmaydp[3] = S1P_obs;
            break;
    }
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
