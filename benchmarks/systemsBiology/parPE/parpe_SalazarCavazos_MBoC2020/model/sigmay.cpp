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
namespace model_SalazarCavazos_MBoC2020 {

void sigmay_SalazarCavazos_MBoC2020(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_pY1068_percent__0_29__dose1 = noiseParameter1_observable_pY1068_percent__0_29__dose1;  // sigmay[0]
    sigma_observable_pY1068_percent__0_66__dose2 = noiseParameter1_observable_pY1068_percent__0_66__dose2;  // sigmay[1]
    sigma_observable_pY1068_percent__0_67__EGF25nM = noiseParameter1_observable_pY1068_percent__0_67__EGF25nM;  // sigmay[2]
    sigma_observable_pY1068_percent__0_69__EGF25nM = noiseParameter1_observable_pY1068_percent__0_69__EGF25nM;  // sigmay[3]
    sigma_observable_pY1068_percent__0_85__EGF25nM = noiseParameter1_observable_pY1068_percent__0_85__EGF25nM;  // sigmay[4]
    sigma_observable_pY1068_percent__0_98__dose3 = noiseParameter1_observable_pY1068_percent__0_98__dose3;  // sigmay[5]
    sigma_observable_pY1173_percent__0_31__dose1 = noiseParameter1_observable_pY1173_percent__0_31__dose1;  // sigmay[6]
    sigma_observable_pY1173_percent__0_66__dose2 = noiseParameter1_observable_pY1173_percent__0_66__dose2;  // sigmay[7]
    sigma_observable_pY1173_percent__0_79__EGF25nM = noiseParameter1_observable_pY1173_percent__0_79__EGF25nM;  // sigmay[8]
    sigma_observable_pY1173_percent__0_8__dose3 = noiseParameter1_observable_pY1173_percent__0_8__dose3;  // sigmay[9]
    sigma_observable_pY1173_percent__0_88__EGF25nM = noiseParameter1_observable_pY1173_percent__0_88__EGF25nM;  // sigmay[10]
    sigma_observable_pY1173_percent__1_02__EGF25nM = noiseParameter1_observable_pY1173_percent__1_02__EGF25nM;  // sigmay[11]
    sigma_observable_phosR_per__0_2__EGF25nM = noiseParameter1_observable_phosR_per__0_2__EGF25nM;  // sigmay[12]
    sigma_observable_phosR_per__0_26__EGF25nM = noiseParameter1_observable_phosR_per__0_26__EGF25nM;  // sigmay[13]
    sigma_observable_phosR_per__0_34__EGF25nM = noiseParameter1_observable_phosR_per__0_34__EGF25nM;  // sigmay[14]
    sigma_observable_phosR_per__0_61__dose1 = noiseParameter1_observable_phosR_per__0_61__dose1;  // sigmay[15]
    sigma_observable_phosR_per__0_97__dose2 = noiseParameter1_observable_phosR_per__0_97__dose2;  // sigmay[16]
    sigma_observable_phosR_per__1_3__dose3 = noiseParameter1_observable_phosR_per__1_3__dose3;  // sigmay[17]
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
