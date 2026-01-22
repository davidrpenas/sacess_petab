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

void Jy_Smith_BMCSystBiol2013(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PI3K_activity__2B, 2)) + 0.5*std::pow(PI3K_activity__2B - mPI3K_activity__2B, 2)/std::pow(sigma_PI3K_activity__2B, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Glucose_uptake__2B, 2)) + 0.5*std::pow(Glucose_uptake__2B - mGlucose_uptake__2B, 2)/std::pow(sigma_Glucose_uptake__2B, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Cell_Bound_Ins__2B, 2)) + 0.5*std::pow(Cell_Bound_Ins__2B - mCell_Bound_Ins__2B, 2)/std::pow(sigma_Cell_Bound_Ins__2B, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_IRSYp__2C, 2)) + 0.5*std::pow(IRSYp__2C - mIRSYp__2C, 2)/std::pow(sigma_IRSYp__2C, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_PTP_activ__2D, 2)) + 0.5*std::pow(PTP_activ__2D - mPTP_activ__2D, 2)/std::pow(sigma_PTP_activ__2D, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Glucose_uptake__120__3B, 2)) + 0.5*std::pow(Glucose_uptake__120__3B - mGlucose_uptake__120__3B, 2)/std::pow(sigma_Glucose_uptake__120__3B, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Glucose_uptake__240__3B, 2)) + 0.5*std::pow(Glucose_uptake__240__3B - mGlucose_uptake__240__3B, 2)/std::pow(sigma_Glucose_uptake__240__3B, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_MnSOD_fold_induction__3C, 2)) + 0.5*std::pow(MnSOD_fold_induction__3C - mMnSOD_fold_induction__3C, 2)/std::pow(sigma_MnSOD_fold_induction__3C, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FOXO4__3C, 2)) + 0.5*std::pow(FOXO4__3C - mFOXO4__3C, 2)/std::pow(sigma_FOXO4__3C, 2);
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
