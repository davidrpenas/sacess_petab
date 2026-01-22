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
namespace model_Okuonghae_ChaosSolitonsFractals2020 {

void Jy_Okuonghae_ChaosSolitonsFractals2020(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Cumulative, 2)) + 0.5*std::pow(-mobservable_Cumulative + observable_Cumulative, 2)/std::pow(sigma_observable_Cumulative, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_ActiveCases, 2)) + 0.5*std::pow(-mobservable_ActiveCases + observable_ActiveCases, 2)/std::pow(sigma_observable_ActiveCases, 2);
            break;
    }
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
