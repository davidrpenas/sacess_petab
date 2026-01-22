#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

static constexpr std::array<std::array<sunindextype, 19>, 18> dJydy_colptrs_SalazarCavazos_MBoC2020_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_SalazarCavazos_MBoC2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_SalazarCavazos_MBoC2020_[index]));
}
} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

static constexpr std::array<std::array<sunindextype, 1>, 18> dJydy_rowvals_SalazarCavazos_MBoC2020_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_SalazarCavazos_MBoC2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_SalazarCavazos_MBoC2020_[index]));
}
} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici




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
#include "dJydy.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void dJydy_SalazarCavazos_MBoC2020(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_pY1068_percent__0_29__dose1 + 1.0*observable_pY1068_percent__0_29__dose1)/std::pow(sigma_observable_pY1068_percent__0_29__dose1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_pY1068_percent__0_66__dose2 + 1.0*observable_pY1068_percent__0_66__dose2)/std::pow(sigma_observable_pY1068_percent__0_66__dose2, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_pY1068_percent__0_67__EGF25nM + 1.0*observable_pY1068_percent__0_67__EGF25nM)/std::pow(sigma_observable_pY1068_percent__0_67__EGF25nM, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_pY1068_percent__0_69__EGF25nM + 1.0*observable_pY1068_percent__0_69__EGF25nM)/std::pow(sigma_observable_pY1068_percent__0_69__EGF25nM, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobservable_pY1068_percent__0_85__EGF25nM + 1.0*observable_pY1068_percent__0_85__EGF25nM)/std::pow(sigma_observable_pY1068_percent__0_85__EGF25nM, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobservable_pY1068_percent__0_98__dose3 + 1.0*observable_pY1068_percent__0_98__dose3)/std::pow(sigma_observable_pY1068_percent__0_98__dose3, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mobservable_pY1173_percent__0_31__dose1 + 1.0*observable_pY1173_percent__0_31__dose1)/std::pow(sigma_observable_pY1173_percent__0_31__dose1, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mobservable_pY1173_percent__0_66__dose2 + 1.0*observable_pY1173_percent__0_66__dose2)/std::pow(sigma_observable_pY1173_percent__0_66__dose2, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mobservable_pY1173_percent__0_79__EGF25nM + 1.0*observable_pY1173_percent__0_79__EGF25nM)/std::pow(sigma_observable_pY1173_percent__0_79__EGF25nM, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*mobservable_pY1173_percent__0_8__dose3 + 1.0*observable_pY1173_percent__0_8__dose3)/std::pow(sigma_observable_pY1173_percent__0_8__dose3, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*mobservable_pY1173_percent__0_88__EGF25nM + 1.0*observable_pY1173_percent__0_88__EGF25nM)/std::pow(sigma_observable_pY1173_percent__0_88__EGF25nM, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*mobservable_pY1173_percent__1_02__EGF25nM + 1.0*observable_pY1173_percent__1_02__EGF25nM)/std::pow(sigma_observable_pY1173_percent__1_02__EGF25nM, 2);
            break;
        case 12:
            dJydy[0] = (-1.0*mobservable_phosR_per__0_2__EGF25nM + 1.0*observable_phosR_per__0_2__EGF25nM)/std::pow(sigma_observable_phosR_per__0_2__EGF25nM, 2);
            break;
        case 13:
            dJydy[0] = (-1.0*mobservable_phosR_per__0_26__EGF25nM + 1.0*observable_phosR_per__0_26__EGF25nM)/std::pow(sigma_observable_phosR_per__0_26__EGF25nM, 2);
            break;
        case 14:
            dJydy[0] = (-1.0*mobservable_phosR_per__0_34__EGF25nM + 1.0*observable_phosR_per__0_34__EGF25nM)/std::pow(sigma_observable_phosR_per__0_34__EGF25nM, 2);
            break;
        case 15:
            dJydy[0] = (-1.0*mobservable_phosR_per__0_61__dose1 + 1.0*observable_phosR_per__0_61__dose1)/std::pow(sigma_observable_phosR_per__0_61__dose1, 2);
            break;
        case 16:
            dJydy[0] = (-1.0*mobservable_phosR_per__0_97__dose2 + 1.0*observable_phosR_per__0_97__dose2)/std::pow(sigma_observable_phosR_per__0_97__dose2, 2);
            break;
        case 17:
            dJydy[0] = (-1.0*mobservable_phosR_per__1_3__dose3 + 1.0*observable_phosR_per__1_3__dose3)/std::pow(sigma_observable_phosR_per__1_3__dose3, 2);
            break;
    }
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
