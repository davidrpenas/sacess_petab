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
namespace model_Bertozzi_PNAS2020 {

void Jy_Bertozzi_PNAS2020(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_y_I_NY, 2)) + 0.5*std::pow(-my_I_NY + y_I_NY, 2)/std::pow(sigma_y_I_NY, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_y_I_CA, 2)) + 0.5*std::pow(-my_I_CA + y_I_CA, 2)/std::pow(sigma_y_I_CA, 2);
            break;
    }
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
