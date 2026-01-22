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

void dJydsigma_Weber_BMC2015(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_yCERTpRN24 - 1.0*std::pow(-myCERTpRN24 + yCERTpRN24, 2)/std::pow(sigma_yCERTpRN24, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_yCERTt - 1.0*std::pow(-myCERTt + yCERTt, 2)/std::pow(sigma_yCERTt, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_yPI4K3BpRN24 - 1.0*std::pow(-myPI4K3BpRN24 + yPI4K3BpRN24, 2)/std::pow(sigma_yPI4K3BpRN24, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_yPI4K3Bt - 1.0*std::pow(-myPI4K3Bt + yPI4K3Bt, 2)/std::pow(sigma_yPI4K3Bt, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_yPKDpN0 - 1.0*std::pow(-myPKDpN0 + yPKDpN0, 2)/std::pow(sigma_yPKDpN0, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_yPKDpN24 - 1.0*std::pow(-myPKDpN24 + yPKDpN24, 2)/std::pow(sigma_yPKDpN24, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_yPKDpN25 - 1.0*std::pow(-myPKDpN25 + yPKDpN25, 2)/std::pow(sigma_yPKDpN25, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_yPKDt - 1.0*std::pow(-myPKDt + yPKDt, 2)/std::pow(sigma_yPKDt, 3);
            break;
    }
}

} // namespace model_Weber_BMC2015
} // namespace amici
