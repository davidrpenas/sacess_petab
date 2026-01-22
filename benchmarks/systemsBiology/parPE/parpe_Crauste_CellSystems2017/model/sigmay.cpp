#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

void sigmay_Crauste_CellSystems2017(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_EarlyEffector__266__model1_data1 = noiseParameter1_observable_EarlyEffector__266__model1_data1;  // sigmay[0]
    sigma_observable_EarlyEffector__526__model1_data1 = noiseParameter1_observable_EarlyEffector__526__model1_data1;  // sigmay[1]
    sigma_observable_EarlyEffector__1540__model1_data1 = noiseParameter1_observable_EarlyEffector__1540__model1_data1;  // sigmay[2]
    sigma_observable_EarlyEffector__2924__model1_data1 = noiseParameter1_observable_EarlyEffector__2924__model1_data1;  // sigmay[3]
    sigma_observable_EarlyEffector__3048__model1_data1 = noiseParameter1_observable_EarlyEffector__3048__model1_data1;  // sigmay[4]
    sigma_observable_EarlyEffector__19722__model1_data1 = noiseParameter1_observable_EarlyEffector__19722__model1_data1;  // sigmay[5]
    sigma_observable_EarlyEffector__23855__model1_data1 = noiseParameter1_observable_EarlyEffector__23855__model1_data1;  // sigmay[6]
    sigma_observable_LateEffector__182__model1_data1 = noiseParameter1_observable_LateEffector__182__model1_data1;  // sigmay[7]
    sigma_observable_LateEffector__1976__model1_data1 = noiseParameter1_observable_LateEffector__1976__model1_data1;  // sigmay[8]
    sigma_observable_LateEffector__2705__model1_data1 = noiseParameter1_observable_LateEffector__2705__model1_data1;  // sigmay[9]
    sigma_observable_LateEffector__4105__model1_data1 = noiseParameter1_observable_LateEffector__4105__model1_data1;  // sigmay[10]
    sigma_observable_LateEffector__6841__model1_data1 = noiseParameter1_observable_LateEffector__6841__model1_data1;  // sigmay[11]
    sigma_observable_LateEffector__7266__model1_data1 = noiseParameter1_observable_LateEffector__7266__model1_data1;  // sigmay[12]
    sigma_observable_LateEffector__23601__model1_data1 = noiseParameter1_observable_LateEffector__23601__model1_data1;  // sigmay[13]
    sigma_observable_Memory__80__model1_data1 = noiseParameter1_observable_Memory__80__model1_data1;  // sigmay[14]
    sigma_observable_Memory__4206__model1_data1 = noiseParameter1_observable_Memory__4206__model1_data1;  // sigmay[15]
    sigma_observable_Memory__5745__model1_data1 = noiseParameter1_observable_Memory__5745__model1_data1;  // sigmay[16]
    sigma_observable_Memory__7909__model1_data1 = noiseParameter1_observable_Memory__7909__model1_data1;  // sigmay[17]
    sigma_observable_Memory__15987__model1_data1 = noiseParameter1_observable_Memory__15987__model1_data1;  // sigmay[18]
    sigma_observable_Naive__64__model1_data1 = noiseParameter1_observable_Naive__64__model1_data1;  // sigmay[19]
    sigma_observable_Naive__290__model1_data1 = noiseParameter1_observable_Naive__290__model1_data1;  // sigmay[20]
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
