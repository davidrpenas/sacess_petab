#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void x_rdata_Isensee_JCB2018(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = AC;
    x_rdata[1] = pAC;
    x_rdata[2] = AC_Fsk;
    x_rdata[3] = pAC_Fsk;
    x_rdata[4] = cAMP;
    x_rdata[5] = Csub;
    x_rdata[6] = Csub_H89;
    x_rdata[7] = RII_C_2;
    x_rdata[8] = RIIp_C_2;
    x_rdata[9] = RIIp_cAMP_C_2;
    x_rdata[10] = RIIp_cAMP_2;
    x_rdata[11] = Sp8_Br_cAMPS;
    x_rdata[12] = Rp8_Br_cAMPS;
    x_rdata[13] = Rp8_pCPT_cAMPS;
    x_rdata[14] = Rp_cAMPS;
    x_rdata[15] = RIIp_Rp8_Br_cAMPS_C_2;
    x_rdata[16] = RIIp_Rp8_pCPT_cAMPS_C_2;
    x_rdata[17] = RIIp_Rp_cAMPS_C_2;
    x_rdata[18] = RIIp_Sp8_Br_cAMPS_C_2;
    x_rdata[19] = RIIp_Sp8_Br_cAMPS_2;
    x_rdata[20] = RIIp_2;
    x_rdata[21] = RII_2;
    x_rdata[22] = IBMX;
    x_rdata[23] = PDE;
    x_rdata[24] = pPDE;
}

} // namespace model_Isensee_JCB2018
} // namespace amici
