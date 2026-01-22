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
namespace model_Perelson_Science1996 {

void Jy_Perelson_Science1996(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mtask0_model0_perelson1_V, 2)*std::pow(sigma_task0_model0_perelson1_V, 2)*std::pow(M_LN10, 2)) + 0.5*std::pow(-std::log(mtask0_model0_perelson1_V)/M_LN10 + std::log(task0_model0_perelson1_V)/M_LN10, 2)/std::pow(sigma_task0_model0_perelson1_V, 2);
            break;
    }
}

} // namespace model_Perelson_Science1996
} // namespace amici
