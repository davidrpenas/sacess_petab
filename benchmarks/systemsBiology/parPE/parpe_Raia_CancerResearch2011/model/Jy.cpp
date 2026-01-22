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

void Jy_Raia_CancerResearch2011(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_CD274mRNA, 2)) + 0.5*std::pow(-mobservable_CD274mRNA + observable_CD274mRNA, 2)/std::pow(sigma_observable_CD274mRNA, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_IL13_cell, 2)) + 0.5*std::pow(-mobservable_IL13_cell + observable_IL13_cell, 2)/std::pow(sigma_observable_IL13_cell, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_RecSurf, 2)) + 0.5*std::pow(-mobservable_RecSurf + observable_RecSurf, 2)/std::pow(sigma_observable_RecSurf, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_SOCS3, 2)) + 0.5*std::pow(-mobservable_SOCS3 + observable_SOCS3, 2)/std::pow(sigma_observable_SOCS3, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_SOCS3mRNA, 2)) + 0.5*std::pow(-mobservable_SOCS3mRNA + observable_SOCS3mRNA, 2)/std::pow(sigma_observable_SOCS3mRNA, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pIL4Ra, 2)) + 0.5*std::pow(-mobservable_pIL4Ra + observable_pIL4Ra, 2)/std::pow(sigma_observable_pIL4Ra, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pJAK2, 2)) + 0.5*std::pow(-mobservable_pJAK2 + observable_pJAK2, 2)/std::pow(sigma_observable_pJAK2, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pSTAT5, 2)) + 0.5*std::pow(-mobservable_pSTAT5 + observable_pSTAT5, 2)/std::pow(sigma_observable_pSTAT5, 2);
            break;
    }
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
