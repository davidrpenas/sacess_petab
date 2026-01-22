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
namespace model_Isensee_JCB2018 {

void Jy_Isensee_JCB2018(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pRII_Microscopy, 2)) + 0.5*std::pow(-mpRII_Microscopy + pRII_Microscopy, 2)/std::pow(sigma_pRII_Microscopy, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_pRII_Western, 2)) + 0.5*std::pow(-mpRII_Western + pRII_Western, 2)/std::pow(sigma_pRII_Western, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_Calpha_Microscopy, 2)) + 0.5*std::pow(Calpha_Microscopy - mCalpha_Microscopy, 2)/std::pow(sigma_Calpha_Microscopy, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_KD_Fsk, 2)) + 0.5*std::pow(-mprior_KD_Fsk + prior_KD_Fsk, 2)/std::pow(sigma_prior_KD_Fsk, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_KD_H89, 2)) + 0.5*std::pow(-mprior_KD_H89 + prior_KD_H89, 2)/std::pow(sigma_prior_KD_H89, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_KD_IBMX, 2)) + 0.5*std::pow(-mprior_KD_IBMX + prior_KD_IBMX, 2)/std::pow(sigma_prior_KD_IBMX, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_KD_cAMP, 2)) + 0.5*std::pow(-mprior_KD_cAMP + prior_KD_cAMP, 2)/std::pow(sigma_prior_KD_cAMP, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_b_Calpha_global, 2)) + 0.5*std::pow(-mprior_b_Calpha_global + prior_b_Calpha_global, 2)/std::pow(sigma_prior_b_Calpha_global, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_b_pRII_global, 2)) + 0.5*std::pow(-mprior_b_pRII_global + prior_b_pRII_global, 2)/std::pow(sigma_prior_b_pRII_global, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_s_pRII_JI09_150302_Drg345_343_CycNuc, 2)) + 0.5*std::pow(-mprior_s_pRII_JI09_150302_Drg345_343_CycNuc + prior_s_pRII_JI09_150302_Drg345_343_CycNuc, 2)/std::pow(sigma_prior_s_pRII_JI09_150302_Drg345_343_CycNuc, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_s_pRII_JI09_150330_Drg350_348_CycNuc, 2)) + 0.5*std::pow(-mprior_s_pRII_JI09_150330_Drg350_348_CycNuc + prior_s_pRII_JI09_150330_Drg350_348_CycNuc, 2)/std::pow(sigma_prior_s_pRII_JI09_150330_Drg350_348_CycNuc, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_s_pRII_JI09_150330_Drg353_351_CycNuc, 2)) + 0.5*std::pow(-mprior_s_pRII_JI09_150330_Drg353_351_CycNuc + prior_s_pRII_JI09_150330_Drg353_351_CycNuc, 2)/std::pow(sigma_prior_s_pRII_JI09_150330_Drg353_351_CycNuc, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_s_pRII_JI09_151102_Drg421_418_Age, 2)) + 0.5*std::pow(-mprior_s_pRII_JI09_151102_Drg421_418_Age + prior_s_pRII_JI09_151102_Drg421_418_Age, 2)/std::pow(sigma_prior_s_pRII_JI09_151102_Drg421_418_Age, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_b_Rp8_Br_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_b_Rp8_Br_cAMPS + prior_xi_b_Rp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_b_Rp8_Br_cAMPS, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_b_Rp8_pCPT_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_b_Rp8_pCPT_cAMPS + prior_xi_b_Rp8_pCPT_cAMPS, 2)/std::pow(sigma_prior_xi_b_Rp8_pCPT_cAMPS, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_b_Rp_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_b_Rp_cAMPS + prior_xi_b_Rp_cAMPS, 2)/std::pow(sigma_prior_xi_b_Rp_cAMPS, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_b_Sp8_Br_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_b_Sp8_Br_cAMPS + prior_xi_b_Sp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_b_Sp8_Br_cAMPS, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kf_RII_2__RII_C_2, 2)) + 0.5*std::pow(-mprior_xi_kf_RII_2__RII_C_2 + prior_xi_kf_RII_2__RII_C_2, 2)/std::pow(sigma_prior_xi_kf_RII_2__RII_C_2, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_kf_RII_C_2__RII_2, 2)) + 0.5*std::pow(-mprior_xi_kf_RII_C_2__RII_2 + prior_xi_kf_RII_C_2__RII_2, 2)/std::pow(sigma_prior_xi_kf_RII_C_2__RII_2, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_KD_Rp8_Br_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_KD_Rp8_Br_cAMPS + prior_xi_KD_Rp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Rp8_Br_cAMPS, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_KD_Rp8_pCPT_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_KD_Rp8_pCPT_cAMPS + prior_xi_KD_Rp8_pCPT_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Rp8_pCPT_cAMPS, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_KD_Rp_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_KD_Rp_cAMPS + prior_xi_KD_Rp_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Rp_cAMPS, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_prior_xi_KD_Sp8_Br_cAMPS, 2)) + 0.5*std::pow(-mprior_xi_KD_Sp8_Br_cAMPS + prior_xi_KD_Sp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Sp8_Br_cAMPS, 2);
            break;
    }
}

} // namespace model_Isensee_JCB2018
} // namespace amici
