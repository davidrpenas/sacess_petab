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
namespace model_Borghans_BiophysChem1997 {

void Jy_Borghans_BiophysChem1997(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mCa, 2)*std::pow(sigma_Ca, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(std::log(Ca)/M_LN10 - std::log(mCa)/M_LN10, 2)/std::pow(sigma_Ca, 2);
            break;
    }
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
