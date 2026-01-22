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
namespace model_Zheng_PNAS2012 {

void Jy_Zheng_PNAS2012(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me0K36me0, 2)) + 0.5*std::pow(-mobservable_K27me0K36me0 + observable_K27me0K36me0, 2)/std::pow(sigma_observable_K27me0K36me0, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me0K36me1, 2)) + 0.5*std::pow(-mobservable_K27me0K36me1 + observable_K27me0K36me1, 2)/std::pow(sigma_observable_K27me0K36me1, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me0K36me2, 2)) + 0.5*std::pow(-mobservable_K27me0K36me2 + observable_K27me0K36me2, 2)/std::pow(sigma_observable_K27me0K36me2, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me0K36me3, 2)) + 0.5*std::pow(-mobservable_K27me0K36me3 + observable_K27me0K36me3, 2)/std::pow(sigma_observable_K27me0K36me3, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me1K36me0, 2)) + 0.5*std::pow(-mobservable_K27me1K36me0 + observable_K27me1K36me0, 2)/std::pow(sigma_observable_K27me1K36me0, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me1K36me1, 2)) + 0.5*std::pow(-mobservable_K27me1K36me1 + observable_K27me1K36me1, 2)/std::pow(sigma_observable_K27me1K36me1, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me1K36me2, 2)) + 0.5*std::pow(-mobservable_K27me1K36me2 + observable_K27me1K36me2, 2)/std::pow(sigma_observable_K27me1K36me2, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me1K36me3, 2)) + 0.5*std::pow(-mobservable_K27me1K36me3 + observable_K27me1K36me3, 2)/std::pow(sigma_observable_K27me1K36me3, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me2K36me0, 2)) + 0.5*std::pow(-mobservable_K27me2K36me0 + observable_K27me2K36me0, 2)/std::pow(sigma_observable_K27me2K36me0, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me2K36me1, 2)) + 0.5*std::pow(-mobservable_K27me2K36me1 + observable_K27me2K36me1, 2)/std::pow(sigma_observable_K27me2K36me1, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me2K36me2, 2)) + 0.5*std::pow(-mobservable_K27me2K36me2 + observable_K27me2K36me2, 2)/std::pow(sigma_observable_K27me2K36me2, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me2K36me3, 2)) + 0.5*std::pow(-mobservable_K27me2K36me3 + observable_K27me2K36me3, 2)/std::pow(sigma_observable_K27me2K36me3, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me3K36me0, 2)) + 0.5*std::pow(-mobservable_K27me3K36me0 + observable_K27me3K36me0, 2)/std::pow(sigma_observable_K27me3K36me0, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me3K36me1, 2)) + 0.5*std::pow(-mobservable_K27me3K36me1 + observable_K27me3K36me1, 2)/std::pow(sigma_observable_K27me3K36me1, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_K27me3K36me2, 2)) + 0.5*std::pow(-mobservable_K27me3K36me2 + observable_K27me3K36me2, 2)/std::pow(sigma_observable_K27me3K36me2, 2);
            break;
    }
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
