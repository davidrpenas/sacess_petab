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
namespace model_Weber_BMC2015 {

void Jy_Weber_BMC2015(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCERTpRN24, 2)) + 0.5*std::pow(-myCERTpRN24 + yCERTpRN24, 2)/std::pow(sigma_yCERTpRN24, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yCERTt, 2)) + 0.5*std::pow(-myCERTt + yCERTt, 2)/std::pow(sigma_yCERTt, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI4K3BpRN24, 2)) + 0.5*std::pow(-myPI4K3BpRN24 + yPI4K3BpRN24, 2)/std::pow(sigma_yPI4K3BpRN24, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPI4K3Bt, 2)) + 0.5*std::pow(-myPI4K3Bt + yPI4K3Bt, 2)/std::pow(sigma_yPI4K3Bt, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPKDpN0, 2)) + 0.5*std::pow(-myPKDpN0 + yPKDpN0, 2)/std::pow(sigma_yPKDpN0, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPKDpN24, 2)) + 0.5*std::pow(-myPKDpN24 + yPKDpN24, 2)/std::pow(sigma_yPKDpN24, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPKDpN25, 2)) + 0.5*std::pow(-myPKDpN25 + yPKDpN25, 2)/std::pow(sigma_yPKDpN25, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yPKDt, 2)) + 0.5*std::pow(-myPKDt + yPKDt, 2)/std::pow(sigma_yPKDt, 2);
            break;
    }
}

} // namespace model_Weber_BMC2015
} // namespace amici
