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
namespace model_Rahman_MBS2016 {

void dJydsigma_Rahman_MBS2016(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_prevalence - 1.0*std::pow(-mobservable_prevalence + observable_prevalence, 2)/std::pow(sigma_observable_prevalence, 3);
            break;
    }
}

} // namespace model_Rahman_MBS2016
} // namespace amici
