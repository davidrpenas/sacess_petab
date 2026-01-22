#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void y_Isensee_JCB2018(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = observableParameter1_pRII_Microscopy*(b_pRII_global + s_pRII_global*(2*RIIp_2 + 2*RIIp_Sp8_Br_cAMPS_2 + 2*RIIp_cAMP_2 + 2*rel_open*(RIIp_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 + RIIp_Rp_cAMPS_C_2) + (2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2)*(rel_open - xi_rel_open*(rel_open - 1))));
    y[1] = s_pRII_Western*(2*RIIp_2 + 2*RIIp_C_2 + 2*RIIp_Rp8_Br_cAMPS_C_2 + 2*RIIp_Rp8_pCPT_cAMPS_C_2 + 2*RIIp_Rp_cAMPS_C_2 + 2*RIIp_Sp8_Br_cAMPS_2 + 2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_2 + 2*RIIp_cAMP_C_2);
    y[2] = observableParameter1_Calpha_Microscopy*(b_Calpha_global + s_Calpha_global*(2*Csub + 2*Csub_H89 + 2*rel_open*(RIIp_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 + RIIp_Rp_cAMPS_C_2) + (2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2)*(rel_open - xi_rel_open*(rel_open - 1))));
    y[3] = std::log(KD_Fsk)/M_LN10;
    y[4] = std::log(KD_H89)/M_LN10;
    y[5] = std::log(KD_IBMX)/M_LN10;
    y[6] = std::log(KD_cAMP)/M_LN10;
    y[7] = std::log(b_Calpha_global)/M_LN10;
    y[8] = std::log(b_pRII_global)/M_LN10;
    y[9] = std::log(s_pRII_JI09_150302_Drg345_343_CycNuc)/M_LN10;
    y[10] = std::log(s_pRII_JI09_150330_Drg350_348_CycNuc)/M_LN10;
    y[11] = std::log(s_pRII_JI09_150330_Drg353_351_CycNuc)/M_LN10;
    y[12] = std::log(s_pRII_JI09_151102_Drg421_418_Age)/M_LN10;
    y[13] = std::log(xi_b_Rp8_Br_cAMPS)/M_LN10;
    y[14] = std::log(xi_b_Rp8_pCPT_cAMPS)/M_LN10;
    y[15] = std::log(xi_b_Rp_cAMPS)/M_LN10;
    y[16] = std::log(xi_b_Sp8_Br_cAMPS)/M_LN10;
    y[17] = std::log(xi_kf_RII_2__RII_C_2)/M_LN10;
    y[18] = std::log(xi_kf_RII_C_2__RII_2)/M_LN10;
    y[19] = std::log(xi_KD_Rp8_Br_cAMPS)/M_LN10;
    y[20] = std::log(xi_KD_Rp8_pCPT_cAMPS)/M_LN10;
    y[21] = std::log(xi_KD_Rp_cAMPS)/M_LN10;
    y[22] = std::log(xi_KD_Sp8_Br_cAMPS)/M_LN10;
}

} // namespace model_Isensee_JCB2018
} // namespace amici
