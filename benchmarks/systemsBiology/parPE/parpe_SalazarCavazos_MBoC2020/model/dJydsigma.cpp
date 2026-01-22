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
namespace model_SalazarCavazos_MBoC2020 {

void dJydsigma_SalazarCavazos_MBoC2020(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_pY1068_percent__0_29__dose1 - 1.0*std::pow(-mobservable_pY1068_percent__0_29__dose1 + observable_pY1068_percent__0_29__dose1, 2)/std::pow(sigma_observable_pY1068_percent__0_29__dose1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_pY1068_percent__0_66__dose2 - 1.0*std::pow(-mobservable_pY1068_percent__0_66__dose2 + observable_pY1068_percent__0_66__dose2, 2)/std::pow(sigma_observable_pY1068_percent__0_66__dose2, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_pY1068_percent__0_67__EGF25nM - 1.0*std::pow(-mobservable_pY1068_percent__0_67__EGF25nM + observable_pY1068_percent__0_67__EGF25nM, 2)/std::pow(sigma_observable_pY1068_percent__0_67__EGF25nM, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_pY1068_percent__0_69__EGF25nM - 1.0*std::pow(-mobservable_pY1068_percent__0_69__EGF25nM + observable_pY1068_percent__0_69__EGF25nM, 2)/std::pow(sigma_observable_pY1068_percent__0_69__EGF25nM, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_observable_pY1068_percent__0_85__EGF25nM - 1.0*std::pow(-mobservable_pY1068_percent__0_85__EGF25nM + observable_pY1068_percent__0_85__EGF25nM, 2)/std::pow(sigma_observable_pY1068_percent__0_85__EGF25nM, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_observable_pY1068_percent__0_98__dose3 - 1.0*std::pow(-mobservable_pY1068_percent__0_98__dose3 + observable_pY1068_percent__0_98__dose3, 2)/std::pow(sigma_observable_pY1068_percent__0_98__dose3, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_observable_pY1173_percent__0_31__dose1 - 1.0*std::pow(-mobservable_pY1173_percent__0_31__dose1 + observable_pY1173_percent__0_31__dose1, 2)/std::pow(sigma_observable_pY1173_percent__0_31__dose1, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_observable_pY1173_percent__0_66__dose2 - 1.0*std::pow(-mobservable_pY1173_percent__0_66__dose2 + observable_pY1173_percent__0_66__dose2, 2)/std::pow(sigma_observable_pY1173_percent__0_66__dose2, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_observable_pY1173_percent__0_79__EGF25nM - 1.0*std::pow(-mobservable_pY1173_percent__0_79__EGF25nM + observable_pY1173_percent__0_79__EGF25nM, 2)/std::pow(sigma_observable_pY1173_percent__0_79__EGF25nM, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_observable_pY1173_percent__0_8__dose3 - 1.0*std::pow(-mobservable_pY1173_percent__0_8__dose3 + observable_pY1173_percent__0_8__dose3, 2)/std::pow(sigma_observable_pY1173_percent__0_8__dose3, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_observable_pY1173_percent__0_88__EGF25nM - 1.0*std::pow(-mobservable_pY1173_percent__0_88__EGF25nM + observable_pY1173_percent__0_88__EGF25nM, 2)/std::pow(sigma_observable_pY1173_percent__0_88__EGF25nM, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_observable_pY1173_percent__1_02__EGF25nM - 1.0*std::pow(-mobservable_pY1173_percent__1_02__EGF25nM + observable_pY1173_percent__1_02__EGF25nM, 2)/std::pow(sigma_observable_pY1173_percent__1_02__EGF25nM, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_observable_phosR_per__0_2__EGF25nM - 1.0*std::pow(-mobservable_phosR_per__0_2__EGF25nM + observable_phosR_per__0_2__EGF25nM, 2)/std::pow(sigma_observable_phosR_per__0_2__EGF25nM, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_observable_phosR_per__0_26__EGF25nM - 1.0*std::pow(-mobservable_phosR_per__0_26__EGF25nM + observable_phosR_per__0_26__EGF25nM, 2)/std::pow(sigma_observable_phosR_per__0_26__EGF25nM, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_observable_phosR_per__0_34__EGF25nM - 1.0*std::pow(-mobservable_phosR_per__0_34__EGF25nM + observable_phosR_per__0_34__EGF25nM, 2)/std::pow(sigma_observable_phosR_per__0_34__EGF25nM, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_observable_phosR_per__0_61__dose1 - 1.0*std::pow(-mobservable_phosR_per__0_61__dose1 + observable_phosR_per__0_61__dose1, 2)/std::pow(sigma_observable_phosR_per__0_61__dose1, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_observable_phosR_per__0_97__dose2 - 1.0*std::pow(-mobservable_phosR_per__0_97__dose2 + observable_phosR_per__0_97__dose2, 2)/std::pow(sigma_observable_phosR_per__0_97__dose2, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_observable_phosR_per__1_3__dose3 - 1.0*std::pow(-mobservable_phosR_per__1_3__dose3 + observable_phosR_per__1_3__dose3, 2)/std::pow(sigma_observable_phosR_per__1_3__dose3, 3);
            break;
    }
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
