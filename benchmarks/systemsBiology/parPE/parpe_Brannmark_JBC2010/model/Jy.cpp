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
namespace model_Brannmark_JBC2010 {

void Jy_Brannmark_JBC2010(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_IR1_P, 2)) + 0.5*std::pow(IR1_P - mIR1_P, 2)/std::pow(sigma_IR1_P, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_IRS1_P, 2)) + 0.5*std::pow(IRS1_P - mIRS1_P, 2)/std::pow(sigma_IRS1_P, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_IRS1_P_DosR, 2)) + 0.5*std::pow(IRS1_P_DosR - mIRS1_P_DosR, 2)/std::pow(sigma_IRS1_P_DosR, 2);
            break;
    }
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
