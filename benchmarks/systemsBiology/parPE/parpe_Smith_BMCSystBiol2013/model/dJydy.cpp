#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

static constexpr std::array<std::array<sunindextype, 10>, 9> dJydy_colptrs_Smith_BMCSystBiol2013_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Smith_BMCSystBiol2013(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Smith_BMCSystBiol2013_[index]));
}
} // namespace model_Smith_BMCSystBiol2013
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

static constexpr std::array<std::array<sunindextype, 1>, 9> dJydy_rowvals_Smith_BMCSystBiol2013_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Smith_BMCSystBiol2013(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Smith_BMCSystBiol2013_[index]));
}
} // namespace model_Smith_BMCSystBiol2013
} // namespace amici




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
#include "dJydy.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void dJydy_Smith_BMCSystBiol2013(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*PI3K_activity__2B - 1.0*mPI3K_activity__2B)/std::pow(sigma_PI3K_activity__2B, 2);
            break;
        case 1:
            dJydy[0] = (1.0*Glucose_uptake__2B - 1.0*mGlucose_uptake__2B)/std::pow(sigma_Glucose_uptake__2B, 2);
            break;
        case 2:
            dJydy[0] = (1.0*Cell_Bound_Ins__2B - 1.0*mCell_Bound_Ins__2B)/std::pow(sigma_Cell_Bound_Ins__2B, 2);
            break;
        case 3:
            dJydy[0] = (1.0*IRSYp__2C - 1.0*mIRSYp__2C)/std::pow(sigma_IRSYp__2C, 2);
            break;
        case 4:
            dJydy[0] = (1.0*PTP_activ__2D - 1.0*mPTP_activ__2D)/std::pow(sigma_PTP_activ__2D, 2);
            break;
        case 5:
            dJydy[0] = (1.0*Glucose_uptake__120__3B - 1.0*mGlucose_uptake__120__3B)/std::pow(sigma_Glucose_uptake__120__3B, 2);
            break;
        case 6:
            dJydy[0] = (1.0*Glucose_uptake__240__3B - 1.0*mGlucose_uptake__240__3B)/std::pow(sigma_Glucose_uptake__240__3B, 2);
            break;
        case 7:
            dJydy[0] = (1.0*MnSOD_fold_induction__3C - 1.0*mMnSOD_fold_induction__3C)/std::pow(sigma_MnSOD_fold_induction__3C, 2);
            break;
        case 8:
            dJydy[0] = (1.0*FOXO4__3C - 1.0*mFOXO4__3C)/std::pow(sigma_FOXO4__3C, 2);
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
