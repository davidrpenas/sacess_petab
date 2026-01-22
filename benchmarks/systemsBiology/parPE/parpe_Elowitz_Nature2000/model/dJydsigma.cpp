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
namespace model_Elowitz_Nature2000 {

void dJydsigma_Elowitz_Nature2000(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_fluorescence - 1.0*std::pow(std::log(fluorescence)/M_LN10 - std::log(mfluorescence)/M_LN10, 2)/std::pow(sigma_fluorescence, 3);
            break;
    }
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
