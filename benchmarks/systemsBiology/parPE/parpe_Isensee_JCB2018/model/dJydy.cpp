#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<std::array<sunindextype, 24>, 23> dJydy_colptrs_Isensee_JCB2018_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Isensee_JCB2018(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Isensee_JCB2018_[index]));
}
} // namespace model_Isensee_JCB2018
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Isensee_JCB2018 {

static constexpr std::array<std::array<sunindextype, 1>, 23> dJydy_rowvals_Isensee_JCB2018_ = {{
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
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Isensee_JCB2018(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Isensee_JCB2018_[index]));
}
} // namespace model_Isensee_JCB2018
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
namespace model_Isensee_JCB2018 {

void dJydy_Isensee_JCB2018(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mpRII_Microscopy + 1.0*pRII_Microscopy)/std::pow(sigma_pRII_Microscopy, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mpRII_Western + 1.0*pRII_Western)/std::pow(sigma_pRII_Western, 2);
            break;
        case 2:
            dJydy[0] = (1.0*Calpha_Microscopy - 1.0*mCalpha_Microscopy)/std::pow(sigma_Calpha_Microscopy, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mprior_KD_Fsk + 1.0*prior_KD_Fsk)/std::pow(sigma_prior_KD_Fsk, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mprior_KD_H89 + 1.0*prior_KD_H89)/std::pow(sigma_prior_KD_H89, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mprior_KD_IBMX + 1.0*prior_KD_IBMX)/std::pow(sigma_prior_KD_IBMX, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mprior_KD_cAMP + 1.0*prior_KD_cAMP)/std::pow(sigma_prior_KD_cAMP, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mprior_b_Calpha_global + 1.0*prior_b_Calpha_global)/std::pow(sigma_prior_b_Calpha_global, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mprior_b_pRII_global + 1.0*prior_b_pRII_global)/std::pow(sigma_prior_b_pRII_global, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*mprior_s_pRII_JI09_150302_Drg345_343_CycNuc + 1.0*prior_s_pRII_JI09_150302_Drg345_343_CycNuc)/std::pow(sigma_prior_s_pRII_JI09_150302_Drg345_343_CycNuc, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*mprior_s_pRII_JI09_150330_Drg350_348_CycNuc + 1.0*prior_s_pRII_JI09_150330_Drg350_348_CycNuc)/std::pow(sigma_prior_s_pRII_JI09_150330_Drg350_348_CycNuc, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*mprior_s_pRII_JI09_150330_Drg353_351_CycNuc + 1.0*prior_s_pRII_JI09_150330_Drg353_351_CycNuc)/std::pow(sigma_prior_s_pRII_JI09_150330_Drg353_351_CycNuc, 2);
            break;
        case 12:
            dJydy[0] = (-1.0*mprior_s_pRII_JI09_151102_Drg421_418_Age + 1.0*prior_s_pRII_JI09_151102_Drg421_418_Age)/std::pow(sigma_prior_s_pRII_JI09_151102_Drg421_418_Age, 2);
            break;
        case 13:
            dJydy[0] = (-1.0*mprior_xi_b_Rp8_Br_cAMPS + 1.0*prior_xi_b_Rp8_Br_cAMPS)/std::pow(sigma_prior_xi_b_Rp8_Br_cAMPS, 2);
            break;
        case 14:
            dJydy[0] = (-1.0*mprior_xi_b_Rp8_pCPT_cAMPS + 1.0*prior_xi_b_Rp8_pCPT_cAMPS)/std::pow(sigma_prior_xi_b_Rp8_pCPT_cAMPS, 2);
            break;
        case 15:
            dJydy[0] = (-1.0*mprior_xi_b_Rp_cAMPS + 1.0*prior_xi_b_Rp_cAMPS)/std::pow(sigma_prior_xi_b_Rp_cAMPS, 2);
            break;
        case 16:
            dJydy[0] = (-1.0*mprior_xi_b_Sp8_Br_cAMPS + 1.0*prior_xi_b_Sp8_Br_cAMPS)/std::pow(sigma_prior_xi_b_Sp8_Br_cAMPS, 2);
            break;
        case 17:
            dJydy[0] = (-1.0*mprior_xi_kf_RII_2__RII_C_2 + 1.0*prior_xi_kf_RII_2__RII_C_2)/std::pow(sigma_prior_xi_kf_RII_2__RII_C_2, 2);
            break;
        case 18:
            dJydy[0] = (-1.0*mprior_xi_kf_RII_C_2__RII_2 + 1.0*prior_xi_kf_RII_C_2__RII_2)/std::pow(sigma_prior_xi_kf_RII_C_2__RII_2, 2);
            break;
        case 19:
            dJydy[0] = (-1.0*mprior_xi_KD_Rp8_Br_cAMPS + 1.0*prior_xi_KD_Rp8_Br_cAMPS)/std::pow(sigma_prior_xi_KD_Rp8_Br_cAMPS, 2);
            break;
        case 20:
            dJydy[0] = (-1.0*mprior_xi_KD_Rp8_pCPT_cAMPS + 1.0*prior_xi_KD_Rp8_pCPT_cAMPS)/std::pow(sigma_prior_xi_KD_Rp8_pCPT_cAMPS, 2);
            break;
        case 21:
            dJydy[0] = (-1.0*mprior_xi_KD_Rp_cAMPS + 1.0*prior_xi_KD_Rp_cAMPS)/std::pow(sigma_prior_xi_KD_Rp_cAMPS, 2);
            break;
        case 22:
            dJydy[0] = (-1.0*mprior_xi_KD_Sp8_Br_cAMPS + 1.0*prior_xi_KD_Sp8_Br_cAMPS)/std::pow(sigma_prior_xi_KD_Sp8_Br_cAMPS, 2);
            break;
    }
}

} // namespace model_Isensee_JCB2018
} // namespace amici
