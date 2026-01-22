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
namespace model_Raia_CancerResearch2011 {

void dJydsigma_Raia_CancerResearch2011(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_CD274mRNA - 1.0*std::pow(-mobservable_CD274mRNA + observable_CD274mRNA, 2)/std::pow(sigma_observable_CD274mRNA, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_IL13_cell - 1.0*std::pow(-mobservable_IL13_cell + observable_IL13_cell, 2)/std::pow(sigma_observable_IL13_cell, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_RecSurf - 1.0*std::pow(-mobservable_RecSurf + observable_RecSurf, 2)/std::pow(sigma_observable_RecSurf, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_SOCS3 - 1.0*std::pow(-mobservable_SOCS3 + observable_SOCS3, 2)/std::pow(sigma_observable_SOCS3, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_SOCS3mRNA - 1.0*std::pow(-mobservable_SOCS3mRNA + observable_SOCS3mRNA, 2)/std::pow(sigma_observable_SOCS3mRNA, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_observable_pIL4Ra - 1.0*std::pow(-mobservable_pIL4Ra + observable_pIL4Ra, 2)/std::pow(sigma_observable_pIL4Ra, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_observable_pJAK2 - 1.0*std::pow(-mobservable_pJAK2 + observable_pJAK2, 2)/std::pow(sigma_observable_pJAK2, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_observable_pSTAT5 - 1.0*std::pow(-mobservable_pSTAT5 + observable_pSTAT5, 2)/std::pow(sigma_observable_pSTAT5, 3);
            break;
    }
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
