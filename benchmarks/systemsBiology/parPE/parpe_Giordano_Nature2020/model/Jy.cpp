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

void Jy_Giordano_Nature2020(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_TotalCases, 2)) + 0.5*std::pow(-mobservable_TotalCases + observable_TotalCases, 2)/std::pow(sigma_observable_TotalCases, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Deaths, 2)) + 0.5*std::pow(-mobservable_Deaths + observable_Deaths, 2)/std::pow(sigma_observable_Deaths, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Recovered, 2)) + 0.5*std::pow(-mobservable_Recovered + observable_Recovered, 2)/std::pow(sigma_observable_Recovered, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_CurrentCases, 2)) + 0.5*std::pow(-mobservable_CurrentCases + observable_CurrentCases, 2)/std::pow(sigma_observable_CurrentCases, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_DiagHome, 2)) + 0.5*std::pow(-mobservable_DiagHome + observable_DiagHome, 2)/std::pow(sigma_observable_DiagHome, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Hospit, 2)) + 0.5*std::pow(-mobservable_Hospit + observable_Hospit, 2)/std::pow(sigma_observable_Hospit, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_ICU, 2)) + 0.5*std::pow(-mobservable_ICU + observable_ICU, 2)/std::pow(sigma_observable_ICU, 2);
            break;
    }
}

} // namespace model_Giordano_Nature2020
} // namespace amici
