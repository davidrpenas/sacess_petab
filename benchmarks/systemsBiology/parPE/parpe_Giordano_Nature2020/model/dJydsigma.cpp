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
namespace model_Giordano_Nature2020 {

void dJydsigma_Giordano_Nature2020(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_TotalCases - 1.0*std::pow(-mobservable_TotalCases + observable_TotalCases, 2)/std::pow(sigma_observable_TotalCases, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_Deaths - 1.0*std::pow(-mobservable_Deaths + observable_Deaths, 2)/std::pow(sigma_observable_Deaths, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_Recovered - 1.0*std::pow(-mobservable_Recovered + observable_Recovered, 2)/std::pow(sigma_observable_Recovered, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_CurrentCases - 1.0*std::pow(-mobservable_CurrentCases + observable_CurrentCases, 2)/std::pow(sigma_observable_CurrentCases, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_DiagHome - 1.0*std::pow(-mobservable_DiagHome + observable_DiagHome, 2)/std::pow(sigma_observable_DiagHome, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_observable_Hospit - 1.0*std::pow(-mobservable_Hospit + observable_Hospit, 2)/std::pow(sigma_observable_Hospit, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_observable_ICU - 1.0*std::pow(-mobservable_ICU + observable_ICU, 2)/std::pow(sigma_observable_ICU, 3);
            break;
    }
}

} // namespace model_Giordano_Nature2020
} // namespace amici
