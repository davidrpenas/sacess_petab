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

void Jy_Crauste_CellSystems2017(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_EarlyEffector__266__model1_data1, 2)) + 0.5*std::pow(-mobservable_EarlyEffector__266__model1_data1 + observable_EarlyEffector__266__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__266__model1_data1, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_EarlyEffector__526__model1_data1, 2)) + 0.5*std::pow(-mobservable_EarlyEffector__526__model1_data1 + observable_EarlyEffector__526__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__526__model1_data1, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_EarlyEffector__1540__model1_data1, 2)) + 0.5*std::pow(-mobservable_EarlyEffector__1540__model1_data1 + observable_EarlyEffector__1540__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__1540__model1_data1, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_EarlyEffector__2924__model1_data1, 2)) + 0.5*std::pow(-mobservable_EarlyEffector__2924__model1_data1 + observable_EarlyEffector__2924__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__2924__model1_data1, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_EarlyEffector__3048__model1_data1, 2)) + 0.5*std::pow(-mobservable_EarlyEffector__3048__model1_data1 + observable_EarlyEffector__3048__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__3048__model1_data1, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_EarlyEffector__19722__model1_data1, 2)) + 0.5*std::pow(-mobservable_EarlyEffector__19722__model1_data1 + observable_EarlyEffector__19722__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__19722__model1_data1, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_EarlyEffector__23855__model1_data1, 2)) + 0.5*std::pow(-mobservable_EarlyEffector__23855__model1_data1 + observable_EarlyEffector__23855__model1_data1, 2)/std::pow(sigma_observable_EarlyEffector__23855__model1_data1, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_LateEffector__182__model1_data1, 2)) + 0.5*std::pow(-mobservable_LateEffector__182__model1_data1 + observable_LateEffector__182__model1_data1, 2)/std::pow(sigma_observable_LateEffector__182__model1_data1, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_LateEffector__1976__model1_data1, 2)) + 0.5*std::pow(-mobservable_LateEffector__1976__model1_data1 + observable_LateEffector__1976__model1_data1, 2)/std::pow(sigma_observable_LateEffector__1976__model1_data1, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_LateEffector__2705__model1_data1, 2)) + 0.5*std::pow(-mobservable_LateEffector__2705__model1_data1 + observable_LateEffector__2705__model1_data1, 2)/std::pow(sigma_observable_LateEffector__2705__model1_data1, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_LateEffector__4105__model1_data1, 2)) + 0.5*std::pow(-mobservable_LateEffector__4105__model1_data1 + observable_LateEffector__4105__model1_data1, 2)/std::pow(sigma_observable_LateEffector__4105__model1_data1, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_LateEffector__6841__model1_data1, 2)) + 0.5*std::pow(-mobservable_LateEffector__6841__model1_data1 + observable_LateEffector__6841__model1_data1, 2)/std::pow(sigma_observable_LateEffector__6841__model1_data1, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_LateEffector__7266__model1_data1, 2)) + 0.5*std::pow(-mobservable_LateEffector__7266__model1_data1 + observable_LateEffector__7266__model1_data1, 2)/std::pow(sigma_observable_LateEffector__7266__model1_data1, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_LateEffector__23601__model1_data1, 2)) + 0.5*std::pow(-mobservable_LateEffector__23601__model1_data1 + observable_LateEffector__23601__model1_data1, 2)/std::pow(sigma_observable_LateEffector__23601__model1_data1, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Memory__80__model1_data1, 2)) + 0.5*std::pow(-mobservable_Memory__80__model1_data1 + observable_Memory__80__model1_data1, 2)/std::pow(sigma_observable_Memory__80__model1_data1, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Memory__4206__model1_data1, 2)) + 0.5*std::pow(-mobservable_Memory__4206__model1_data1 + observable_Memory__4206__model1_data1, 2)/std::pow(sigma_observable_Memory__4206__model1_data1, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Memory__5745__model1_data1, 2)) + 0.5*std::pow(-mobservable_Memory__5745__model1_data1 + observable_Memory__5745__model1_data1, 2)/std::pow(sigma_observable_Memory__5745__model1_data1, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Memory__7909__model1_data1, 2)) + 0.5*std::pow(-mobservable_Memory__7909__model1_data1 + observable_Memory__7909__model1_data1, 2)/std::pow(sigma_observable_Memory__7909__model1_data1, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Memory__15987__model1_data1, 2)) + 0.5*std::pow(-mobservable_Memory__15987__model1_data1 + observable_Memory__15987__model1_data1, 2)/std::pow(sigma_observable_Memory__15987__model1_data1, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Naive__64__model1_data1, 2)) + 0.5*std::pow(-mobservable_Naive__64__model1_data1 + observable_Naive__64__model1_data1, 2)/std::pow(sigma_observable_Naive__64__model1_data1, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_Naive__290__model1_data1, 2)) + 0.5*std::pow(-mobservable_Naive__290__model1_data1 + observable_Naive__290__model1_data1, 2)/std::pow(sigma_observable_Naive__290__model1_data1, 2);
            break;
    }
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
