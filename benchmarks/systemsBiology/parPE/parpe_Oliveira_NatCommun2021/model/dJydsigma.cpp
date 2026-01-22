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
namespace model_Oliveira_NatCommun2021 {

void dJydsigma_Oliveira_NatCommun2021(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_cumulative_deaths - 1.0*std::pow(cumulative_deaths - mcumulative_deaths, 2)/std::pow(sigma_cumulative_deaths, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_cumulative_cases - 1.0*std::pow(cumulative_cases - mcumulative_cases, 2)/std::pow(sigma_cumulative_cases, 3);
            break;
    }
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
