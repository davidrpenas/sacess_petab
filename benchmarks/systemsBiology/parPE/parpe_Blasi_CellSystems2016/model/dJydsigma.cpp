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
namespace model_Blasi_CellSystems2016 {

void dJydsigma_Blasi_CellSystems2016(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_0ac - 1.0*std::pow(-std::log(mobservable_0ac) + std::log(observable_0ac), 2)/std::pow(sigma_observable_0ac, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_4ac - 1.0*std::pow(-std::log(mobservable_4ac) + std::log(observable_4ac), 2)/std::pow(sigma_observable_4ac, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_k12 - 1.0*std::pow(-std::log(mobservable_k12) + std::log(observable_k12), 2)/std::pow(sigma_observable_k12, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_k12k16 - 1.0*std::pow(-std::log(mobservable_k12k16) + std::log(observable_k12k16), 2)/std::pow(sigma_observable_k12k16, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_k16 - 1.0*std::pow(-std::log(mobservable_k16) + std::log(observable_k16), 2)/std::pow(sigma_observable_k16, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_observable_k5 - 1.0*std::pow(-std::log(mobservable_k5) + std::log(observable_k5), 2)/std::pow(sigma_observable_k5, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_observable_k5k12 - 1.0*std::pow(-std::log(mobservable_k5k12) + std::log(observable_k5k12), 2)/std::pow(sigma_observable_k5k12, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_observable_k5k12k16 - 1.0*std::pow(-std::log(mobservable_k5k12k16) + std::log(observable_k5k12k16), 2)/std::pow(sigma_observable_k5k12k16, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_observable_k5k16 - 1.0*std::pow(-std::log(mobservable_k5k16) + std::log(observable_k5k16), 2)/std::pow(sigma_observable_k5k16, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_observable_k5k8 - 1.0*std::pow(-std::log(mobservable_k5k8) + std::log(observable_k5k8), 2)/std::pow(sigma_observable_k5k8, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_observable_k5k8k12 - 1.0*std::pow(-std::log(mobservable_k5k8k12) + std::log(observable_k5k8k12), 2)/std::pow(sigma_observable_k5k8k12, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_observable_k5k8k16 - 1.0*std::pow(-std::log(mobservable_k5k8k16) + std::log(observable_k5k8k16), 2)/std::pow(sigma_observable_k5k8k16, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_observable_k8 - 1.0*std::pow(-std::log(mobservable_k8) + std::log(observable_k8), 2)/std::pow(sigma_observable_k8, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_observable_k8k12 - 1.0*std::pow(-std::log(mobservable_k8k12) + std::log(observable_k8k12), 2)/std::pow(sigma_observable_k8k12, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_observable_k8k12k16 - 1.0*std::pow(-std::log(mobservable_k8k12k16) + std::log(observable_k8k12k16), 2)/std::pow(sigma_observable_k8k12k16, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_observable_k8k16 - 1.0*std::pow(-std::log(mobservable_k8k16) + std::log(observable_k8k16), 2)/std::pow(sigma_observable_k8k16, 3);
            break;
    }
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
