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
namespace model_Smith_BMCSystBiol2013 {

void dJydsigma_Smith_BMCSystBiol2013(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_PI3K_activity__2B - 1.0*std::pow(PI3K_activity__2B - mPI3K_activity__2B, 2)/std::pow(sigma_PI3K_activity__2B, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_Glucose_uptake__2B - 1.0*std::pow(Glucose_uptake__2B - mGlucose_uptake__2B, 2)/std::pow(sigma_Glucose_uptake__2B, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_Cell_Bound_Ins__2B - 1.0*std::pow(Cell_Bound_Ins__2B - mCell_Bound_Ins__2B, 2)/std::pow(sigma_Cell_Bound_Ins__2B, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_IRSYp__2C - 1.0*std::pow(IRSYp__2C - mIRSYp__2C, 2)/std::pow(sigma_IRSYp__2C, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_PTP_activ__2D - 1.0*std::pow(PTP_activ__2D - mPTP_activ__2D, 2)/std::pow(sigma_PTP_activ__2D, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_Glucose_uptake__120__3B - 1.0*std::pow(Glucose_uptake__120__3B - mGlucose_uptake__120__3B, 2)/std::pow(sigma_Glucose_uptake__120__3B, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_Glucose_uptake__240__3B - 1.0*std::pow(Glucose_uptake__240__3B - mGlucose_uptake__240__3B, 2)/std::pow(sigma_Glucose_uptake__240__3B, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_MnSOD_fold_induction__3C - 1.0*std::pow(MnSOD_fold_induction__3C - mMnSOD_fold_induction__3C, 2)/std::pow(sigma_MnSOD_fold_induction__3C, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_FOXO4__3C - 1.0*std::pow(FOXO4__3C - mFOXO4__3C, 2)/std::pow(sigma_FOXO4__3C, 3);
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
