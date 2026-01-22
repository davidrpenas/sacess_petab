#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<std::array<sunindextype, 5>, 4> dJydy_colptrs_Armistead_CellDeathDis2024_ = {{
    {0, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1}, 
    {0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Armistead_CellDeathDis2024(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Armistead_CellDeathDis2024_[index]));
}
} // namespace model_Armistead_CellDeathDis2024
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<std::array<sunindextype, 1>, 4> dJydy_rowvals_Armistead_CellDeathDis2024_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Armistead_CellDeathDis2024(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Armistead_CellDeathDis2024_[index]));
}
} // namespace model_Armistead_CellDeathDis2024
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
namespace model_Armistead_CellDeathDis2024 {

void dJydy_Armistead_CellDeathDis2024(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*Sphinga_obs - 1.0*mSphinga_obs)/std::pow(sigma_Sphinga_obs, 2);
            break;
        case 1:
            dJydy[0] = (1.0*Cer_obs - 1.0*mCer_obs)/std::pow(sigma_Cer_obs, 2);
            break;
        case 2:
            dJydy[0] = (1.0*Sphingo_obs - 1.0*mSphingo_obs)/std::pow(sigma_Sphingo_obs, 2);
            break;
        case 3:
            dJydy[0] = (1.0*S1P_obs - 1.0*mS1P_obs)/std::pow(sigma_S1P_obs, 2);
            break;
    }
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
