#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Blasi_CellSystems2016 {

void sigmay_Blasi_CellSystems2016(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_0ac = noiseParameter1_observable_0ac;  // sigmay[0]
    sigma_observable_4ac = noiseParameter1_observable_4ac;  // sigmay[1]
    sigma_observable_k12 = noiseParameter1_observable_k12;  // sigmay[2]
    sigma_observable_k12k16 = noiseParameter1_observable_k12k16;  // sigmay[3]
    sigma_observable_k16 = noiseParameter1_observable_k16;  // sigmay[4]
    sigma_observable_k5 = noiseParameter1_observable_k5;  // sigmay[5]
    sigma_observable_k5k12 = noiseParameter1_observable_k5k12;  // sigmay[6]
    sigma_observable_k5k12k16 = noiseParameter1_observable_k5k12k16;  // sigmay[7]
    sigma_observable_k5k16 = noiseParameter1_observable_k5k16;  // sigmay[8]
    sigma_observable_k5k8 = noiseParameter1_observable_k5k8;  // sigmay[9]
    sigma_observable_k5k8k12 = noiseParameter1_observable_k5k8k12;  // sigmay[10]
    sigma_observable_k5k8k16 = noiseParameter1_observable_k5k8k16;  // sigmay[11]
    sigma_observable_k8 = noiseParameter1_observable_k8;  // sigmay[12]
    sigma_observable_k8k12 = noiseParameter1_observable_k8k12;  // sigmay[13]
    sigma_observable_k8k12k16 = noiseParameter1_observable_k8k12k16;  // sigmay[14]
    sigma_observable_k8k16 = noiseParameter1_observable_k8k16;  // sigmay[15]
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
