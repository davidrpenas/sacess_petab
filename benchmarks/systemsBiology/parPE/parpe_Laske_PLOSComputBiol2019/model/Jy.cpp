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

void Jy_Laske_PLOSComputBiol2019(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mRM5, 2)*std::pow(sigma_RM5, 2)) + 0.5*std::pow(std::log(RM5) - std::log(mRM5), 2)/std::pow(sigma_RM5, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mRVSegTot, 2)*std::pow(sigma_RVSegTot, 2)) + 0.5*std::pow(std::log(RVSegTot) - std::log(mRVSegTot), 2)/std::pow(sigma_RVSegTot, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mRCSegTot, 2)*std::pow(sigma_RCSegTot, 2)) + 0.5*std::pow(std::log(RCSegTot) - std::log(mRCSegTot), 2)/std::pow(sigma_RCSegTot, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mVrel, 2)*std::pow(sigma_Vrel, 2)) + 0.5*std::pow(std::log(Vrel) - std::log(mVrel), 2)/std::pow(sigma_Vrel, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_IntNucOffset, 2)) + 0.5*std::pow(IntNucOffset - mIntNucOffset, 2)/std::pow(sigma_IntNucOffset, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_1, 2)) + 0.5*std::pow(FracNucInt_1 - mFracNucInt_1, 2)/std::pow(sigma_FracNucInt_1, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_2, 2)) + 0.5*std::pow(FracNucInt_2 - mFracNucInt_2, 2)/std::pow(sigma_FracNucInt_2, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_3, 2)) + 0.5*std::pow(FracNucInt_3 - mFracNucInt_3, 2)/std::pow(sigma_FracNucInt_3, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_4, 2)) + 0.5*std::pow(FracNucInt_4 - mFracNucInt_4, 2)/std::pow(sigma_FracNucInt_4, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_5, 2)) + 0.5*std::pow(FracNucInt_5 - mFracNucInt_5, 2)/std::pow(sigma_FracNucInt_5, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_6, 2)) + 0.5*std::pow(FracNucInt_6 - mFracNucInt_6, 2)/std::pow(sigma_FracNucInt_6, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_7, 2)) + 0.5*std::pow(FracNucInt_7 - mFracNucInt_7, 2)/std::pow(sigma_FracNucInt_7, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_FracNucInt_8, 2)) + 0.5*std::pow(FracNucInt_8 - mFracNucInt_8, 2)/std::pow(sigma_FracNucInt_8, 2);
            break;
    }
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
