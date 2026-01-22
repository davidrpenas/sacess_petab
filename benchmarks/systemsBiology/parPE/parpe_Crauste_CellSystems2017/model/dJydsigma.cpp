#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void dJydsigma_Crauste_CellSystems2017(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_EarlyEffector__266__model1_data1 - 1.0*std::pow(-mobservable_EarlyEffector__266__model1_data1 + observable_EarlyEffector__266__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__266__model1_data1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_EarlyEffector__526__model1_data1 - 1.0*std::pow(-mobservable_EarlyEffector__526__model1_data1 + observable_EarlyEffector__526__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__526__model1_data1, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_EarlyEffector__1540__model1_data1 - 1.0*std::pow(-mobservable_EarlyEffector__1540__model1_data1 + observable_EarlyEffector__1540__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__1540__model1_data1, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_EarlyEffector__2924__model1_data1 - 1.0*std::pow(-mobservable_EarlyEffector__2924__model1_data1 + observable_EarlyEffector__2924__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__2924__model1_data1, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_EarlyEffector__3048__model1_data1 - 1.0*std::pow(-mobservable_EarlyEffector__3048__model1_data1 + observable_EarlyEffector__3048__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__3048__model1_data1, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_observable_EarlyEffector__19722__model1_data1 - 1.0*std::pow(-mobservable_EarlyEffector__19722__model1_data1 + observable_EarlyEffector__19722__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__19722__model1_data1, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_observable_EarlyEffector__23855__model1_data1 - 1.0*std::pow(-mobservable_EarlyEffector__23855__model1_data1 + observable_EarlyEffector__23855__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__23855__model1_data1, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_observable_LateEffector__182__model1_data1 - 1.0*std::pow(-mobservable_LateEffector__182__model1_data1 + observable_LateEffector__182__model1_data1, 2)/std::pow(sigma_observable_LateEffector__182__model1_data1, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_observable_LateEffector__1976__model1_data1 - 1.0*std::pow(-mobservable_LateEffector__1976__model1_data1 + observable_LateEffector__1976__model1_data1, 2)/std::pow(sigma_observable_LateEffector__1976__model1_data1, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_observable_LateEffector__2705__model1_data1 - 1.0*std::pow(-mobservable_LateEffector__2705__model1_data1 + observable_LateEffector__2705__model1_data1, 2)/std::pow(sigma_observable_LateEffector__2705__model1_data1, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_observable_LateEffector__4105__model1_data1 - 1.0*std::pow(-mobservable_LateEffector__4105__model1_data1 + observable_LateEffector__4105__model1_data1, 2)/std::pow(sigma_observable_LateEffector__4105__model1_data1, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_observable_LateEffector__6841__model1_data1 - 1.0*std::pow(-mobservable_LateEffector__6841__model1_data1 + observable_LateEffector__6841__model1_data1, 2)/std::pow(sigma_observable_LateEffector__6841__model1_data1, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_observable_LateEffector__7266__model1_data1 - 1.0*std::pow(-mobservable_LateEffector__7266__model1_data1 + observable_LateEffector__7266__model1_data1, 2)/std::pow(sigma_observable_LateEffector__7266__model1_data1, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_observable_LateEffector__23601__model1_data1 - 1.0*std::pow(-mobservable_LateEffector__23601__model1_data1 + observable_LateEffector__23601__model1_data1, 2)/std::pow(sigma_observable_LateEffector__23601__model1_data1, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_observable_Memory__80__model1_data1 - 1.0*std::pow(-mobservable_Memory__80__model1_data1 + observable_Memory__80__model1_data1, 2)/std::pow(sigma_observable_Memory__80__model1_data1, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_observable_Memory__4206__model1_data1 - 1.0*std::pow(-mobservable_Memory__4206__model1_data1 + observable_Memory__4206__model1_data1, 2)/std::pow(sigma_observable_Memory__4206__model1_data1, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_observable_Memory__5745__model1_data1 - 1.0*std::pow(-mobservable_Memory__5745__model1_data1 + observable_Memory__5745__model1_data1, 2)/std::pow(sigma_observable_Memory__5745__model1_data1, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_observable_Memory__7909__model1_data1 - 1.0*std::pow(-mobservable_Memory__7909__model1_data1 + observable_Memory__7909__model1_data1, 2)/std::pow(sigma_observable_Memory__7909__model1_data1, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_observable_Memory__15987__model1_data1 - 1.0*std::pow(-mobservable_Memory__15987__model1_data1 + observable_Memory__15987__model1_data1, 2)/std::pow(sigma_observable_Memory__15987__model1_data1, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_observable_Naive__64__model1_data1 - 1.0*std::pow(-mobservable_Naive__64__model1_data1 + observable_Naive__64__model1_data1, 2)/std::pow(sigma_observable_Naive__64__model1_data1, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_observable_Naive__290__model1_data1 - 1.0*std::pow(-mobservable_Naive__290__model1_data1 + observable_Naive__290__model1_data1, 2)/std::pow(sigma_observable_Naive__290__model1_data1, 3);
            break;
    }
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
