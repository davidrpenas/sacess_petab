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

void dJydsigma_Perelson_Science1996(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_task0_model0_perelson1_V - 1.0*std::pow(-std::log(mtask0_model0_perelson1_V)/M_LN10 + std::log(task0_model0_perelson1_V)/M_LN10, 2)/std::pow(sigma_task0_model0_perelson1_V, 3);
            break;
    }
}

} // namespace model_Perelson_Science1996
} // namespace amici
