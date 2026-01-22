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

void Jy_Blasi_CellSystems2016(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_0ac, 2)*std::pow(sigma_observable_0ac, 2)) + 0.5*std::pow(-std::log(mobservable_0ac) + std::log(observable_0ac), 2)/std::pow(sigma_observable_0ac, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_4ac, 2)*std::pow(sigma_observable_4ac, 2)) + 0.5*std::pow(-std::log(mobservable_4ac) + std::log(observable_4ac), 2)/std::pow(sigma_observable_4ac, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k12, 2)*std::pow(sigma_observable_k12, 2)) + 0.5*std::pow(-std::log(mobservable_k12) + std::log(observable_k12), 2)/std::pow(sigma_observable_k12, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k12k16, 2)*std::pow(sigma_observable_k12k16, 2)) + 0.5*std::pow(-std::log(mobservable_k12k16) + std::log(observable_k12k16), 2)/std::pow(sigma_observable_k12k16, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k16, 2)*std::pow(sigma_observable_k16, 2)) + 0.5*std::pow(-std::log(mobservable_k16) + std::log(observable_k16), 2)/std::pow(sigma_observable_k16, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k5, 2)*std::pow(sigma_observable_k5, 2)) + 0.5*std::pow(-std::log(mobservable_k5) + std::log(observable_k5), 2)/std::pow(sigma_observable_k5, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k5k12, 2)*std::pow(sigma_observable_k5k12, 2)) + 0.5*std::pow(-std::log(mobservable_k5k12) + std::log(observable_k5k12), 2)/std::pow(sigma_observable_k5k12, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k5k12k16, 2)*std::pow(sigma_observable_k5k12k16, 2)) + 0.5*std::pow(-std::log(mobservable_k5k12k16) + std::log(observable_k5k12k16), 2)/std::pow(sigma_observable_k5k12k16, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k5k16, 2)*std::pow(sigma_observable_k5k16, 2)) + 0.5*std::pow(-std::log(mobservable_k5k16) + std::log(observable_k5k16), 2)/std::pow(sigma_observable_k5k16, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k5k8, 2)*std::pow(sigma_observable_k5k8, 2)) + 0.5*std::pow(-std::log(mobservable_k5k8) + std::log(observable_k5k8), 2)/std::pow(sigma_observable_k5k8, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k5k8k12, 2)*std::pow(sigma_observable_k5k8k12, 2)) + 0.5*std::pow(-std::log(mobservable_k5k8k12) + std::log(observable_k5k8k12), 2)/std::pow(sigma_observable_k5k8k12, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k5k8k16, 2)*std::pow(sigma_observable_k5k8k16, 2)) + 0.5*std::pow(-std::log(mobservable_k5k8k16) + std::log(observable_k5k8k16), 2)/std::pow(sigma_observable_k5k8k16, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k8, 2)*std::pow(sigma_observable_k8, 2)) + 0.5*std::pow(-std::log(mobservable_k8) + std::log(observable_k8), 2)/std::pow(sigma_observable_k8, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k8k12, 2)*std::pow(sigma_observable_k8k12, 2)) + 0.5*std::pow(-std::log(mobservable_k8k12) + std::log(observable_k8k12), 2)/std::pow(sigma_observable_k8k12, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k8k12k16, 2)*std::pow(sigma_observable_k8k12k16, 2)) + 0.5*std::pow(-std::log(mobservable_k8k12k16) + std::log(observable_k8k12k16), 2)/std::pow(sigma_observable_k8k12k16, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(mobservable_k8k16, 2)*std::pow(sigma_observable_k8k16, 2)) + 0.5*std::pow(-std::log(mobservable_k8k16) + std::log(observable_k8k16), 2)/std::pow(sigma_observable_k8k16, 2);
            break;
    }
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
