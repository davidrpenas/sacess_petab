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
namespace model_Zheng_PNAS2012 {

void sigmay_Zheng_PNAS2012(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_K27me0K36me0 = noiseParameter1_observable_K27me0K36me0;  // sigmay[0]
    sigma_observable_K27me0K36me1 = noiseParameter1_observable_K27me0K36me1;  // sigmay[1]
    sigma_observable_K27me0K36me2 = noiseParameter1_observable_K27me0K36me2;  // sigmay[2]
    sigma_observable_K27me0K36me3 = noiseParameter1_observable_K27me0K36me3;  // sigmay[3]
    sigma_observable_K27me1K36me0 = noiseParameter1_observable_K27me1K36me0;  // sigmay[4]
    sigma_observable_K27me1K36me1 = noiseParameter1_observable_K27me1K36me1;  // sigmay[5]
    sigma_observable_K27me1K36me2 = noiseParameter1_observable_K27me1K36me2;  // sigmay[6]
    sigma_observable_K27me1K36me3 = noiseParameter1_observable_K27me1K36me3;  // sigmay[7]
    sigma_observable_K27me2K36me0 = noiseParameter1_observable_K27me2K36me0;  // sigmay[8]
    sigma_observable_K27me2K36me1 = noiseParameter1_observable_K27me2K36me1;  // sigmay[9]
    sigma_observable_K27me2K36me2 = noiseParameter1_observable_K27me2K36me2;  // sigmay[10]
    sigma_observable_K27me2K36me3 = noiseParameter1_observable_K27me2K36me3;  // sigmay[11]
    sigma_observable_K27me3K36me0 = noiseParameter1_observable_K27me3K36me0;  // sigmay[12]
    sigma_observable_K27me3K36me1 = noiseParameter1_observable_K27me3K36me1;  // sigmay[13]
    sigma_observable_K27me3K36me2 = noiseParameter1_observable_K27me3K36me2;  // sigmay[14]
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
