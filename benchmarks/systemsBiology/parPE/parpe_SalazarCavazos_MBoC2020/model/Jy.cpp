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

void Jy_SalazarCavazos_MBoC2020(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1068_percent__0_29__dose1, 2)) + 0.5*std::pow(-mobservable_pY1068_percent__0_29__dose1 + observable_pY1068_percent__0_29__dose1, 2)/std::pow(sigma_observable_pY1068_percent__0_29__dose1, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1068_percent__0_66__dose2, 2)) + 0.5*std::pow(-mobservable_pY1068_percent__0_66__dose2 + observable_pY1068_percent__0_66__dose2, 2)/std::pow(sigma_observable_pY1068_percent__0_66__dose2, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1068_percent__0_67__EGF25nM, 2)) + 0.5*std::pow(-mobservable_pY1068_percent__0_67__EGF25nM + observable_pY1068_percent__0_67__EGF25nM, 2)/std::pow(sigma_observable_pY1068_percent__0_67__EGF25nM, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1068_percent__0_69__EGF25nM, 2)) + 0.5*std::pow(-mobservable_pY1068_percent__0_69__EGF25nM + observable_pY1068_percent__0_69__EGF25nM, 2)/std::pow(sigma_observable_pY1068_percent__0_69__EGF25nM, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1068_percent__0_85__EGF25nM, 2)) + 0.5*std::pow(-mobservable_pY1068_percent__0_85__EGF25nM + observable_pY1068_percent__0_85__EGF25nM, 2)/std::pow(sigma_observable_pY1068_percent__0_85__EGF25nM, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1068_percent__0_98__dose3, 2)) + 0.5*std::pow(-mobservable_pY1068_percent__0_98__dose3 + observable_pY1068_percent__0_98__dose3, 2)/std::pow(sigma_observable_pY1068_percent__0_98__dose3, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1173_percent__0_31__dose1, 2)) + 0.5*std::pow(-mobservable_pY1173_percent__0_31__dose1 + observable_pY1173_percent__0_31__dose1, 2)/std::pow(sigma_observable_pY1173_percent__0_31__dose1, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1173_percent__0_66__dose2, 2)) + 0.5*std::pow(-mobservable_pY1173_percent__0_66__dose2 + observable_pY1173_percent__0_66__dose2, 2)/std::pow(sigma_observable_pY1173_percent__0_66__dose2, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1173_percent__0_79__EGF25nM, 2)) + 0.5*std::pow(-mobservable_pY1173_percent__0_79__EGF25nM + observable_pY1173_percent__0_79__EGF25nM, 2)/std::pow(sigma_observable_pY1173_percent__0_79__EGF25nM, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1173_percent__0_8__dose3, 2)) + 0.5*std::pow(-mobservable_pY1173_percent__0_8__dose3 + observable_pY1173_percent__0_8__dose3, 2)/std::pow(sigma_observable_pY1173_percent__0_8__dose3, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1173_percent__0_88__EGF25nM, 2)) + 0.5*std::pow(-mobservable_pY1173_percent__0_88__EGF25nM + observable_pY1173_percent__0_88__EGF25nM, 2)/std::pow(sigma_observable_pY1173_percent__0_88__EGF25nM, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_pY1173_percent__1_02__EGF25nM, 2)) + 0.5*std::pow(-mobservable_pY1173_percent__1_02__EGF25nM + observable_pY1173_percent__1_02__EGF25nM, 2)/std::pow(sigma_observable_pY1173_percent__1_02__EGF25nM, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_phosR_per__0_2__EGF25nM, 2)) + 0.5*std::pow(-mobservable_phosR_per__0_2__EGF25nM + observable_phosR_per__0_2__EGF25nM, 2)/std::pow(sigma_observable_phosR_per__0_2__EGF25nM, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_phosR_per__0_26__EGF25nM, 2)) + 0.5*std::pow(-mobservable_phosR_per__0_26__EGF25nM + observable_phosR_per__0_26__EGF25nM, 2)/std::pow(sigma_observable_phosR_per__0_26__EGF25nM, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_phosR_per__0_34__EGF25nM, 2)) + 0.5*std::pow(-mobservable_phosR_per__0_34__EGF25nM + observable_phosR_per__0_34__EGF25nM, 2)/std::pow(sigma_observable_phosR_per__0_34__EGF25nM, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_phosR_per__0_61__dose1, 2)) + 0.5*std::pow(-mobservable_phosR_per__0_61__dose1 + observable_phosR_per__0_61__dose1, 2)/std::pow(sigma_observable_phosR_per__0_61__dose1, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_phosR_per__0_97__dose2, 2)) + 0.5*std::pow(-mobservable_phosR_per__0_97__dose2 + observable_phosR_per__0_97__dose2, 2)/std::pow(sigma_observable_phosR_per__0_97__dose2, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_observable_phosR_per__1_3__dose3, 2)) + 0.5*std::pow(-mobservable_phosR_per__1_3__dose3 + observable_phosR_per__1_3__dose3, 2)/std::pow(sigma_observable_phosR_per__1_3__dose3, 2);
            break;
    }
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
