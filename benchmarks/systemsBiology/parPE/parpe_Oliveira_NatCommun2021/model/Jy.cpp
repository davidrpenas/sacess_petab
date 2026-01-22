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

void Jy_Oliveira_NatCommun2021(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_cumulative_deaths, 2)) + 0.5*std::pow(cumulative_deaths - mcumulative_deaths, 2)/std::pow(sigma_cumulative_deaths, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_cumulative_cases, 2)) + 0.5*std::pow(cumulative_cases - mcumulative_cases, 2)/std::pow(sigma_cumulative_cases, 2);
            break;
    }
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
