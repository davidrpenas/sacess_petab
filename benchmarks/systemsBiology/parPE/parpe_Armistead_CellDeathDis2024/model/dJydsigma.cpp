#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

void dJydsigma_Armistead_CellDeathDis2024(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_Sphinga_obs - 1.0*std::pow(Sphinga_obs - mSphinga_obs, 2)/std::pow(sigma_Sphinga_obs, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_Cer_obs - 1.0*std::pow(Cer_obs - mCer_obs, 2)/std::pow(sigma_Cer_obs, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_Sphingo_obs - 1.0*std::pow(Sphingo_obs - mSphingo_obs, 2)/std::pow(sigma_Sphingo_obs, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_S1P_obs - 1.0*std::pow(S1P_obs - mS1P_obs, 2)/std::pow(sigma_S1P_obs, 3);
            break;
    }
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
