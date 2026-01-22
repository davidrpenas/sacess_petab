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
namespace model_Isensee_JCB2018 {

void sigmay_Isensee_JCB2018(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_pRII_Microscopy = noiseParameter1_pRII_Microscopy*rho_pRII_Microscopy*s_pRII_global;  // sigmay[0]
    sigma_pRII_Western = rho_pRII_Western*s_pRII_Western;  // sigmay[1]
    sigma_Calpha_Microscopy = rho_Calpha_Microscopy*s_Calpha_global;  // sigmay[2]
    sigma_prior_KD_Fsk = noiseParameter1_prior_KD_Fsk;  // sigmay[3]
    sigma_prior_KD_H89 = noiseParameter1_prior_KD_H89;  // sigmay[4]
    sigma_prior_KD_IBMX = noiseParameter1_prior_KD_IBMX;  // sigmay[5]
    sigma_prior_KD_cAMP = noiseParameter1_prior_KD_cAMP;  // sigmay[6]
    sigma_prior_b_Calpha_global = noiseParameter1_prior_b_Calpha_global;  // sigmay[7]
    sigma_prior_b_pRII_global = noiseParameter1_prior_b_pRII_global;  // sigmay[8]
    sigma_prior_s_pRII_JI09_150302_Drg345_343_CycNuc = noiseParameter1_prior_s_pRII_JI09_150302_Drg345_343_CycNuc;  // sigmay[9]
    sigma_prior_s_pRII_JI09_150330_Drg350_348_CycNuc = noiseParameter1_prior_s_pRII_JI09_150330_Drg350_348_CycNuc;  // sigmay[10]
    sigma_prior_s_pRII_JI09_150330_Drg353_351_CycNuc = noiseParameter1_prior_s_pRII_JI09_150330_Drg353_351_CycNuc;  // sigmay[11]
    sigma_prior_s_pRII_JI09_151102_Drg421_418_Age = noiseParameter1_prior_s_pRII_JI09_151102_Drg421_418_Age;  // sigmay[12]
    sigma_prior_xi_b_Rp8_Br_cAMPS = noiseParameter1_prior_xi_b_Rp8_Br_cAMPS;  // sigmay[13]
    sigma_prior_xi_b_Rp8_pCPT_cAMPS = noiseParameter1_prior_xi_b_Rp8_pCPT_cAMPS;  // sigmay[14]
    sigma_prior_xi_b_Rp_cAMPS = noiseParameter1_prior_xi_b_Rp_cAMPS;  // sigmay[15]
    sigma_prior_xi_b_Sp8_Br_cAMPS = noiseParameter1_prior_xi_b_Sp8_Br_cAMPS;  // sigmay[16]
    sigma_prior_xi_kf_RII_2__RII_C_2 = noiseParameter1_prior_xi_kf_RII_2__RII_C_2;  // sigmay[17]
    sigma_prior_xi_kf_RII_C_2__RII_2 = noiseParameter1_prior_xi_kf_RII_C_2__RII_2;  // sigmay[18]
    sigma_prior_xi_KD_Rp8_Br_cAMPS = noiseParameter1_prior_xi_KD_Rp8_Br_cAMPS;  // sigmay[19]
    sigma_prior_xi_KD_Rp8_pCPT_cAMPS = noiseParameter1_prior_xi_KD_Rp8_pCPT_cAMPS;  // sigmay[20]
    sigma_prior_xi_KD_Rp_cAMPS = noiseParameter1_prior_xi_KD_Rp_cAMPS;  // sigmay[21]
    sigma_prior_xi_KD_Sp8_Br_cAMPS = noiseParameter1_prior_xi_KD_Sp8_Br_cAMPS;  // sigmay[22]
}

} // namespace model_Isensee_JCB2018
} // namespace amici
