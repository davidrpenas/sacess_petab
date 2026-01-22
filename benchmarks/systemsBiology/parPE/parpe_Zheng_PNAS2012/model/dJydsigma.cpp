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

void dJydsigma_Zheng_PNAS2012(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_K27me0K36me0 - 1.0*std::pow(-mobservable_K27me0K36me0 + observable_K27me0K36me0, 2)/std::pow(sigma_observable_K27me0K36me0, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_K27me0K36me1 - 1.0*std::pow(-mobservable_K27me0K36me1 + observable_K27me0K36me1, 2)/std::pow(sigma_observable_K27me0K36me1, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_K27me0K36me2 - 1.0*std::pow(-mobservable_K27me0K36me2 + observable_K27me0K36me2, 2)/std::pow(sigma_observable_K27me0K36me2, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_K27me0K36me3 - 1.0*std::pow(-mobservable_K27me0K36me3 + observable_K27me0K36me3, 2)/std::pow(sigma_observable_K27me0K36me3, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_K27me1K36me0 - 1.0*std::pow(-mobservable_K27me1K36me0 + observable_K27me1K36me0, 2)/std::pow(sigma_observable_K27me1K36me0, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_observable_K27me1K36me1 - 1.0*std::pow(-mobservable_K27me1K36me1 + observable_K27me1K36me1, 2)/std::pow(sigma_observable_K27me1K36me1, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_observable_K27me1K36me2 - 1.0*std::pow(-mobservable_K27me1K36me2 + observable_K27me1K36me2, 2)/std::pow(sigma_observable_K27me1K36me2, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_observable_K27me1K36me3 - 1.0*std::pow(-mobservable_K27me1K36me3 + observable_K27me1K36me3, 2)/std::pow(sigma_observable_K27me1K36me3, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_observable_K27me2K36me0 - 1.0*std::pow(-mobservable_K27me2K36me0 + observable_K27me2K36me0, 2)/std::pow(sigma_observable_K27me2K36me0, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_observable_K27me2K36me1 - 1.0*std::pow(-mobservable_K27me2K36me1 + observable_K27me2K36me1, 2)/std::pow(sigma_observable_K27me2K36me1, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_observable_K27me2K36me2 - 1.0*std::pow(-mobservable_K27me2K36me2 + observable_K27me2K36me2, 2)/std::pow(sigma_observable_K27me2K36me2, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_observable_K27me2K36me3 - 1.0*std::pow(-mobservable_K27me2K36me3 + observable_K27me2K36me3, 2)/std::pow(sigma_observable_K27me2K36me3, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_observable_K27me3K36me0 - 1.0*std::pow(-mobservable_K27me3K36me0 + observable_K27me3K36me0, 2)/std::pow(sigma_observable_K27me3K36me0, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_observable_K27me3K36me1 - 1.0*std::pow(-mobservable_K27me3K36me1 + observable_K27me3K36me1, 2)/std::pow(sigma_observable_K27me3K36me1, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_observable_K27me3K36me2 - 1.0*std::pow(-mobservable_K27me3K36me2 + observable_K27me3K36me2, 2)/std::pow(sigma_observable_K27me3K36me2, 3);
            break;
    }
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
