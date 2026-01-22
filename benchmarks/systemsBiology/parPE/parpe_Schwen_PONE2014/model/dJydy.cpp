#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<std::array<sunindextype, 11>, 10> dJydy_colptrs_Schwen_PONE2014_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Schwen_PONE2014(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Schwen_PONE2014_[index]));
}
} // namespace model_Schwen_PONE2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Schwen_PONE2014 {

static constexpr std::array<std::array<sunindextype, 1>, 10> dJydy_rowvals_Schwen_PONE2014_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Schwen_PONE2014(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Schwen_PONE2014_[index]));
}
} // namespace model_Schwen_PONE2014
} // namespace amici




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
#include "dJydy.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void dJydy_Schwen_PONE2014(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*std::log(mobservable_IR1)/M_LN10 + 1.0*std::log(observable_IR1)/M_LN10)/(observable_IR1*std::pow(sigma_observable_IR1, 2)*M_LN10);
            break;
        case 1:
            dJydy[0] = (-1.0*std::log(mobservable_IR2)/M_LN10 + 1.0*std::log(observable_IR2)/M_LN10)/(observable_IR2*std::pow(sigma_observable_IR2, 2)*M_LN10);
            break;
        case 2:
            dJydy[0] = (-1.0*std::log(mobservable_IRsum)/M_LN10 + 1.0*std::log(observable_IRsum)/M_LN10)/(observable_IRsum*std::pow(sigma_observable_IRsum, 2)*M_LN10);
            break;
        case 3:
            dJydy[0] = (-1.0*std::log(mobservable_Insulin)/M_LN10 + 1.0*std::log(observable_Insulin)/M_LN10)/(observable_Insulin*std::pow(sigma_observable_Insulin, 2)*M_LN10);
            break;
        case 4:
            dJydy[0] = (-1.0*mprior_ini_R1 + 1.0*prior_ini_R1)/std::pow(sigma_prior_ini_R1, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mprior_ini_R2fold + 1.0*prior_ini_R2fold)/std::pow(sigma_prior_ini_R2fold, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mprior_ka1 + 1.0*prior_ka1)/std::pow(sigma_prior_ka1, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mprior_ka2fold + 1.0*prior_ka2fold)/std::pow(sigma_prior_ka2fold, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mprior_kd1 + 1.0*prior_kd1)/std::pow(sigma_prior_kd1, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*mprior_kd2fold + 1.0*prior_kd2fold)/std::pow(sigma_prior_kd2fold, 2);
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
