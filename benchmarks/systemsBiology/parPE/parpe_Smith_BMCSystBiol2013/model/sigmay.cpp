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
namespace model_Smith_BMCSystBiol2013 {

void sigmay_Smith_BMCSystBiol2013(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_PI3K_activity__2B = 1.0;  // sigmay[0]
    sigma_Glucose_uptake__2B = 1.0;  // sigmay[1]
    sigma_Cell_Bound_Ins__2B = 1.0;  // sigmay[2]
    sigma_IRSYp__2C = 1.0;  // sigmay[3]
    sigma_PTP_activ__2D = 1.0;  // sigmay[4]
    sigma_Glucose_uptake__120__3B = 1.0;  // sigmay[5]
    sigma_Glucose_uptake__240__3B = 1.0;  // sigmay[6]
    sigma_MnSOD_fold_induction__3C = 1.0;  // sigmay[7]
    sigma_FOXO4__3C = 1.0;  // sigmay[8]
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
