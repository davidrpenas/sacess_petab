#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void dJydsigma_Beer_MolBioSystems2014(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_Bacnorm - 1.0*std::pow(Bacnorm - mBacnorm, 2)/std::pow(sigma_Bacnorm, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_IndconcNormRange - 1.0*std::pow(IndconcNormRange - mIndconcNormRange, 2)/std::pow(sigma_IndconcNormRange, 3);
            break;
    }
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
