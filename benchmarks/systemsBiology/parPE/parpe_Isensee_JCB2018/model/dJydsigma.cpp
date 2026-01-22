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

void dJydsigma_Isensee_JCB2018(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_pRII_Microscopy - 1.0*std::pow(-mpRII_Microscopy + pRII_Microscopy, 2)/std::pow(sigma_pRII_Microscopy, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_pRII_Western - 1.0*std::pow(-mpRII_Western + pRII_Western, 2)/std::pow(sigma_pRII_Western, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_Calpha_Microscopy - 1.0*std::pow(Calpha_Microscopy - mCalpha_Microscopy, 2)/std::pow(sigma_Calpha_Microscopy, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_prior_KD_Fsk - 1.0*std::pow(-mprior_KD_Fsk + prior_KD_Fsk, 2)/std::pow(sigma_prior_KD_Fsk, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_prior_KD_H89 - 1.0*std::pow(-mprior_KD_H89 + prior_KD_H89, 2)/std::pow(sigma_prior_KD_H89, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_prior_KD_IBMX - 1.0*std::pow(-mprior_KD_IBMX + prior_KD_IBMX, 2)/std::pow(sigma_prior_KD_IBMX, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_prior_KD_cAMP - 1.0*std::pow(-mprior_KD_cAMP + prior_KD_cAMP, 2)/std::pow(sigma_prior_KD_cAMP, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_prior_b_Calpha_global - 1.0*std::pow(-mprior_b_Calpha_global + prior_b_Calpha_global, 2)/std::pow(sigma_prior_b_Calpha_global, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_prior_b_pRII_global - 1.0*std::pow(-mprior_b_pRII_global + prior_b_pRII_global, 2)/std::pow(sigma_prior_b_pRII_global, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_prior_s_pRII_JI09_150302_Drg345_343_CycNuc - 1.0*std::pow(-mprior_s_pRII_JI09_150302_Drg345_343_CycNuc + prior_s_pRII_JI09_150302_Drg345_343_CycNuc, 2)/std::pow(sigma_prior_s_pRII_JI09_150302_Drg345_343_CycNuc, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_prior_s_pRII_JI09_150330_Drg350_348_CycNuc - 1.0*std::pow(-mprior_s_pRII_JI09_150330_Drg350_348_CycNuc + prior_s_pRII_JI09_150330_Drg350_348_CycNuc, 2)/std::pow(sigma_prior_s_pRII_JI09_150330_Drg350_348_CycNuc, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_prior_s_pRII_JI09_150330_Drg353_351_CycNuc - 1.0*std::pow(-mprior_s_pRII_JI09_150330_Drg353_351_CycNuc + prior_s_pRII_JI09_150330_Drg353_351_CycNuc, 2)/std::pow(sigma_prior_s_pRII_JI09_150330_Drg353_351_CycNuc, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_prior_s_pRII_JI09_151102_Drg421_418_Age - 1.0*std::pow(-mprior_s_pRII_JI09_151102_Drg421_418_Age + prior_s_pRII_JI09_151102_Drg421_418_Age, 2)/std::pow(sigma_prior_s_pRII_JI09_151102_Drg421_418_Age, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_prior_xi_b_Rp8_Br_cAMPS - 1.0*std::pow(-mprior_xi_b_Rp8_Br_cAMPS + prior_xi_b_Rp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_b_Rp8_Br_cAMPS, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_prior_xi_b_Rp8_pCPT_cAMPS - 1.0*std::pow(-mprior_xi_b_Rp8_pCPT_cAMPS + prior_xi_b_Rp8_pCPT_cAMPS, 2)/std::pow(sigma_prior_xi_b_Rp8_pCPT_cAMPS, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_prior_xi_b_Rp_cAMPS - 1.0*std::pow(-mprior_xi_b_Rp_cAMPS + prior_xi_b_Rp_cAMPS, 2)/std::pow(sigma_prior_xi_b_Rp_cAMPS, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_prior_xi_b_Sp8_Br_cAMPS - 1.0*std::pow(-mprior_xi_b_Sp8_Br_cAMPS + prior_xi_b_Sp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_b_Sp8_Br_cAMPS, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_prior_xi_kf_RII_2__RII_C_2 - 1.0*std::pow(-mprior_xi_kf_RII_2__RII_C_2 + prior_xi_kf_RII_2__RII_C_2, 2)/std::pow(sigma_prior_xi_kf_RII_2__RII_C_2, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_prior_xi_kf_RII_C_2__RII_2 - 1.0*std::pow(-mprior_xi_kf_RII_C_2__RII_2 + prior_xi_kf_RII_C_2__RII_2, 2)/std::pow(sigma_prior_xi_kf_RII_C_2__RII_2, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_prior_xi_KD_Rp8_Br_cAMPS - 1.0*std::pow(-mprior_xi_KD_Rp8_Br_cAMPS + prior_xi_KD_Rp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Rp8_Br_cAMPS, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_prior_xi_KD_Rp8_pCPT_cAMPS - 1.0*std::pow(-mprior_xi_KD_Rp8_pCPT_cAMPS + prior_xi_KD_Rp8_pCPT_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Rp8_pCPT_cAMPS, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_prior_xi_KD_Rp_cAMPS - 1.0*std::pow(-mprior_xi_KD_Rp_cAMPS + prior_xi_KD_Rp_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Rp_cAMPS, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_prior_xi_KD_Sp8_Br_cAMPS - 1.0*std::pow(-mprior_xi_KD_Sp8_Br_cAMPS + prior_xi_KD_Sp8_Br_cAMPS, 2)/std::pow(sigma_prior_xi_KD_Sp8_Br_cAMPS, 3);
            break;
    }
}

} // namespace model_Isensee_JCB2018
} // namespace amici
