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
namespace model_Zhao_QuantBiol2020 {

void dJydsigma_Zhao_QuantBiol2020(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_confirmed_infected - 1.0*std::pow(-mobservable_confirmed_infected + observable_confirmed_infected, 2)/std::pow(sigma_observable_confirmed_infected, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_susceptible - 1.0*std::pow(-mobservable_susceptible + observable_susceptible, 2)/std::pow(sigma_observable_susceptible, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_unquarantined_infected - 1.0*std::pow(-mobservable_unquarantined_infected + observable_unquarantined_infected, 2)/std::pow(sigma_observable_unquarantined_infected, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_quarantined_infected - 1.0*std::pow(-mobservable_quarantined_infected + observable_quarantined_infected, 2)/std::pow(sigma_observable_quarantined_infected, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_cumulative_infected - 1.0*std::pow(-mobservable_cumulative_infected + observable_cumulative_infected, 2)/std::pow(sigma_observable_cumulative_infected, 3);
            break;
    }
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
