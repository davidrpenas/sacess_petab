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
namespace model_Laske_PLOSComputBiol2019 {

void sigmay_Laske_PLOSComputBiol2019(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_RM5 = noiseParameter1_RM5;  // sigmay[0]
    sigma_RVSegTot = noiseParameter1_RVSegTot;  // sigmay[1]
    sigma_RCSegTot = noiseParameter1_RCSegTot;  // sigmay[2]
    sigma_Vrel = noiseParameter1_Vrel;  // sigmay[3]
    sigma_IntNucOffset = 13.0334748500416*noiseParameter1_IntNucOffset;  // sigmay[4]
    sigma_FracNucInt_1 = 3.1095605262909198*noiseParameter1_FracNucInt_1;  // sigmay[5]
    sigma_FracNucInt_2 = 3.4465091034262501*noiseParameter1_FracNucInt_2;  // sigmay[6]
    sigma_FracNucInt_3 = 1.9654070825149701*noiseParameter1_FracNucInt_3;  // sigmay[7]
    sigma_FracNucInt_4 = 3.0367677443843699*noiseParameter1_FracNucInt_4;  // sigmay[8]
    sigma_FracNucInt_5 = 4.0031196168321896*noiseParameter1_FracNucInt_5;  // sigmay[9]
    sigma_FracNucInt_6 = 4.8783911282307004*noiseParameter1_FracNucInt_6;  // sigmay[10]
    sigma_FracNucInt_7 = 3.9635379818204601*noiseParameter1_FracNucInt_7;  // sigmay[11]
    sigma_FracNucInt_8 = 4.8820897165046002*noiseParameter1_FracNucInt_8;  // sigmay[12]
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
