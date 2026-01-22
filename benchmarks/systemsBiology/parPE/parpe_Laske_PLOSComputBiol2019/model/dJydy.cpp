#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<std::array<sunindextype, 14>, 13> dJydy_colptrs_Laske_PLOSComputBiol2019_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Laske_PLOSComputBiol2019_[index]));
}
} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<std::array<sunindextype, 1>, 13> dJydy_rowvals_Laske_PLOSComputBiol2019_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
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

void dJydy_rowvals_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Laske_PLOSComputBiol2019_[index]));
}
} // namespace model_Laske_PLOSComputBiol2019
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
namespace model_Laske_PLOSComputBiol2019 {

void dJydy_Laske_PLOSComputBiol2019(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*std::log(RM5) - 1.0*std::log(mRM5))/(RM5*std::pow(sigma_RM5, 2));
            break;
        case 1:
            dJydy[0] = (1.0*std::log(RVSegTot) - 1.0*std::log(mRVSegTot))/(RVSegTot*std::pow(sigma_RVSegTot, 2));
            break;
        case 2:
            dJydy[0] = (1.0*std::log(RCSegTot) - 1.0*std::log(mRCSegTot))/(RCSegTot*std::pow(sigma_RCSegTot, 2));
            break;
        case 3:
            dJydy[0] = (1.0*std::log(Vrel) - 1.0*std::log(mVrel))/(Vrel*std::pow(sigma_Vrel, 2));
            break;
        case 4:
            dJydy[0] = (1.0*IntNucOffset - 1.0*mIntNucOffset)/std::pow(sigma_IntNucOffset, 2);
            break;
        case 5:
            dJydy[0] = (1.0*FracNucInt_1 - 1.0*mFracNucInt_1)/std::pow(sigma_FracNucInt_1, 2);
            break;
        case 6:
            dJydy[0] = (1.0*FracNucInt_2 - 1.0*mFracNucInt_2)/std::pow(sigma_FracNucInt_2, 2);
            break;
        case 7:
            dJydy[0] = (1.0*FracNucInt_3 - 1.0*mFracNucInt_3)/std::pow(sigma_FracNucInt_3, 2);
            break;
        case 8:
            dJydy[0] = (1.0*FracNucInt_4 - 1.0*mFracNucInt_4)/std::pow(sigma_FracNucInt_4, 2);
            break;
        case 9:
            dJydy[0] = (1.0*FracNucInt_5 - 1.0*mFracNucInt_5)/std::pow(sigma_FracNucInt_5, 2);
            break;
        case 10:
            dJydy[0] = (1.0*FracNucInt_6 - 1.0*mFracNucInt_6)/std::pow(sigma_FracNucInt_6, 2);
            break;
        case 11:
            dJydy[0] = (1.0*FracNucInt_7 - 1.0*mFracNucInt_7)/std::pow(sigma_FracNucInt_7, 2);
            break;
        case 12:
            dJydy[0] = (1.0*FracNucInt_8 - 1.0*mFracNucInt_8)/std::pow(sigma_FracNucInt_8, 2);
            break;
    }
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
