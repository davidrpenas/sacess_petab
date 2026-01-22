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
namespace model_Laske_PLOSComputBiol2019 {

void dJydsigma_Laske_PLOSComputBiol2019(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_RM5 - 1.0*std::pow(std::log(RM5) - std::log(mRM5), 2)/std::pow(sigma_RM5, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_RVSegTot - 1.0*std::pow(std::log(RVSegTot) - std::log(mRVSegTot), 2)/std::pow(sigma_RVSegTot, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_RCSegTot - 1.0*std::pow(std::log(RCSegTot) - std::log(mRCSegTot), 2)/std::pow(sigma_RCSegTot, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_Vrel - 1.0*std::pow(std::log(Vrel) - std::log(mVrel), 2)/std::pow(sigma_Vrel, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_IntNucOffset - 1.0*std::pow(IntNucOffset - mIntNucOffset, 2)/std::pow(sigma_IntNucOffset, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_FracNucInt_1 - 1.0*std::pow(FracNucInt_1 - mFracNucInt_1, 2)/std::pow(sigma_FracNucInt_1, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_FracNucInt_2 - 1.0*std::pow(FracNucInt_2 - mFracNucInt_2, 2)/std::pow(sigma_FracNucInt_2, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_FracNucInt_3 - 1.0*std::pow(FracNucInt_3 - mFracNucInt_3, 2)/std::pow(sigma_FracNucInt_3, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_FracNucInt_4 - 1.0*std::pow(FracNucInt_4 - mFracNucInt_4, 2)/std::pow(sigma_FracNucInt_4, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_FracNucInt_5 - 1.0*std::pow(FracNucInt_5 - mFracNucInt_5, 2)/std::pow(sigma_FracNucInt_5, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_FracNucInt_6 - 1.0*std::pow(FracNucInt_6 - mFracNucInt_6, 2)/std::pow(sigma_FracNucInt_6, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_FracNucInt_7 - 1.0*std::pow(FracNucInt_7 - mFracNucInt_7, 2)/std::pow(sigma_FracNucInt_7, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_FracNucInt_8 - 1.0*std::pow(FracNucInt_8 - mFracNucInt_8, 2)/std::pow(sigma_FracNucInt_8, 3);
            break;
    }
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
