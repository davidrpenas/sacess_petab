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
namespace model_Schwen_PONE2014 {

void dJydsigma_Schwen_PONE2014(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_IR1 - 1.0*std::pow(-std::log(mobservable_IR1)/M_LN10 + std::log(observable_IR1)/M_LN10, 2)/std::pow(sigma_observable_IR1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_IR2 - 1.0*std::pow(-std::log(mobservable_IR2)/M_LN10 + std::log(observable_IR2)/M_LN10, 2)/std::pow(sigma_observable_IR2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_IRsum - 1.0*std::pow(-std::log(mobservable_IRsum)/M_LN10 + std::log(observable_IRsum)/M_LN10, 2)/std::pow(sigma_observable_IRsum, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_Insulin - 1.0*std::pow(-std::log(mobservable_Insulin)/M_LN10 + std::log(observable_Insulin)/M_LN10, 2)/std::pow(sigma_observable_Insulin, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_prior_ini_R1 - 1.0*std::pow(-mprior_ini_R1 + prior_ini_R1, 2)/std::pow(sigma_prior_ini_R1, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_prior_ini_R2fold - 1.0*std::pow(-mprior_ini_R2fold + prior_ini_R2fold, 2)/std::pow(sigma_prior_ini_R2fold, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_prior_ka1 - 1.0*std::pow(-mprior_ka1 + prior_ka1, 2)/std::pow(sigma_prior_ka1, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_prior_ka2fold - 1.0*std::pow(-mprior_ka2fold + prior_ka2fold, 2)/std::pow(sigma_prior_ka2fold, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_prior_kd1 - 1.0*std::pow(-mprior_kd1 + prior_kd1, 2)/std::pow(sigma_prior_kd1, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_prior_kd2fold - 1.0*std::pow(-mprior_kd2fold + prior_kd2fold, 2)/std::pow(sigma_prior_kd2fold, 3);
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
