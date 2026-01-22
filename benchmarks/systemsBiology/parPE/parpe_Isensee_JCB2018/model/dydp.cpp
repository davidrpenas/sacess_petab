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

void dydp_Isensee_JCB2018(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 0:
            dydp[3] = 1/(KD_Fsk*M_LN10);
            break;
        case 1:
            dydp[4] = 1/(KD_H89*M_LN10);
            break;
        case 2:
            dydp[5] = 1/(KD_IBMX*M_LN10);
            break;
        case 3:
            dydp[6] = 1/(KD_cAMP*M_LN10);
            break;
        case 22:
            dydp[13] = 1/(xi_b_Rp8_Br_cAMPS*M_LN10);
            break;
        case 23:
            dydp[14] = 1/(xi_b_Rp8_pCPT_cAMPS*M_LN10);
            break;
        case 24:
            dydp[15] = 1/(xi_b_Rp_cAMPS*M_LN10);
            break;
        case 25:
            dydp[16] = 1/(xi_b_Sp8_Br_cAMPS*M_LN10);
            break;
        case 26:
            dydp[17] = 1/(xi_kf_RII_2__RII_C_2*M_LN10);
            break;
        case 27:
            dydp[18] = 1/(xi_kf_RII_C_2__RII_2*M_LN10);
            break;
        case 28:
            dydp[19] = 1/(xi_KD_Rp8_Br_cAMPS*M_LN10);
            break;
        case 29:
            dydp[20] = 1/(xi_KD_Rp8_pCPT_cAMPS*M_LN10);
            break;
        case 30:
            dydp[21] = 1/(xi_KD_Rp_cAMPS*M_LN10);
            break;
        case 31:
            dydp[22] = 1/(xi_KD_Sp8_Br_cAMPS*M_LN10);
            break;
        case 32:
            dydp[0] = observableParameter1_pRII_Microscopy;
            dydp[8] = 1/(b_pRII_global*M_LN10);
            break;
        case 33:
            dydp[0] = b_pRII_global + s_pRII_global*(2*RIIp_2 + 2*RIIp_Sp8_Br_cAMPS_2 + 2*RIIp_cAMP_2 + 2*rel_open*(RIIp_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 + RIIp_Rp_cAMPS_C_2) + (2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2)*(rel_open - xi_rel_open*(rel_open - 1)));
            break;
        case 34:
            dydp[0] = observableParameter1_pRII_Microscopy*s_pRII_global*(2*RIIp_C_2 + 2*RIIp_Rp8_Br_cAMPS_C_2 + 2*RIIp_Rp8_pCPT_cAMPS_C_2 + 2*RIIp_Rp_cAMPS_C_2 + (1 - xi_rel_open)*(2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2));
            dydp[2] = observableParameter1_Calpha_Microscopy*s_Calpha_global*(2*RIIp_C_2 + 2*RIIp_Rp8_Br_cAMPS_C_2 + 2*RIIp_Rp8_pCPT_cAMPS_C_2 + 2*RIIp_Rp_cAMPS_C_2 + (1 - xi_rel_open)*(2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2));
            break;
        case 35:
            dydp[0] = observableParameter1_pRII_Microscopy*(2*RIIp_2 + 2*RIIp_Sp8_Br_cAMPS_2 + 2*RIIp_cAMP_2 + 2*rel_open*(RIIp_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 + RIIp_Rp_cAMPS_C_2) + (2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2)*(rel_open - xi_rel_open*(rel_open - 1)));
            break;
        case 36:
            dydp[0] = observableParameter1_pRII_Microscopy*s_pRII_global*(1 - rel_open)*(2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2);
            dydp[2] = observableParameter1_Calpha_Microscopy*s_Calpha_global*(1 - rel_open)*(2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2);
            break;
        case 37:
            dydp[1] = 2*RIIp_2 + 2*RIIp_C_2 + 2*RIIp_Rp8_Br_cAMPS_C_2 + 2*RIIp_Rp8_pCPT_cAMPS_C_2 + 2*RIIp_Rp_cAMPS_C_2 + 2*RIIp_Sp8_Br_cAMPS_2 + 2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_2 + 2*RIIp_cAMP_C_2;
            break;
        case 38:
            dydp[2] = observableParameter1_Calpha_Microscopy;
            dydp[7] = 1/(b_Calpha_global*M_LN10);
            break;
        case 39:
            dydp[2] = b_Calpha_global + s_Calpha_global*(2*Csub + 2*Csub_H89 + 2*rel_open*(RIIp_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 + RIIp_Rp_cAMPS_C_2) + (2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2)*(rel_open - xi_rel_open*(rel_open - 1)));
            break;
        case 40:
            dydp[2] = observableParameter1_Calpha_Microscopy*(2*Csub + 2*Csub_H89 + 2*rel_open*(RIIp_C_2 + RIIp_Rp8_Br_cAMPS_C_2 + RIIp_Rp8_pCPT_cAMPS_C_2 + RIIp_Rp_cAMPS_C_2) + (2*RIIp_Sp8_Br_cAMPS_C_2 + 2*RIIp_cAMP_C_2)*(rel_open - xi_rel_open*(rel_open - 1)));
            break;
        case 41:
            dydp[9] = 1/(s_pRII_JI09_150302_Drg345_343_CycNuc*M_LN10);
            break;
        case 42:
            dydp[10] = 1/(s_pRII_JI09_150330_Drg350_348_CycNuc*M_LN10);
            break;
        case 43:
            dydp[11] = 1/(s_pRII_JI09_150330_Drg353_351_CycNuc*M_LN10);
            break;
        case 44:
            dydp[12] = 1/(s_pRII_JI09_151102_Drg421_418_Age*M_LN10);
            break;
    }
}

} // namespace model_Isensee_JCB2018
} // namespace amici
