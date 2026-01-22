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

void dJydsigma_Okuonghae_ChaosSolitonsFractals2020(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_Cumulative - 1.0*std::pow(-mobservable_Cumulative + observable_Cumulative, 2)/std::pow(sigma_observable_Cumulative, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_ActiveCases - 1.0*std::pow(-mobservable_ActiveCases + observable_ActiveCases, 2)/std::pow(sigma_observable_ActiveCases, 3);
            break;
    }
}

} // namespace model_Okuonghae_ChaosSolitonsFractals2020
} // namespace amici
