#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Isensee_JCB2018 {

void x_solver_Isensee_JCB2018(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = AC;
    x_solver[1] = pAC;
    x_solver[2] = AC_Fsk;
    x_solver[3] = pAC_Fsk;
    x_solver[4] = cAMP;
    x_solver[5] = Csub;
    x_solver[6] = Csub_H89;
    x_solver[7] = RII_C_2;
    x_solver[8] = RIIp_C_2;
    x_solver[9] = RIIp_cAMP_C_2;
    x_solver[10] = RIIp_cAMP_2;
    x_solver[11] = Sp8_Br_cAMPS;
    x_solver[12] = Rp8_Br_cAMPS;
    x_solver[13] = Rp8_pCPT_cAMPS;
    x_solver[14] = Rp_cAMPS;
    x_solver[15] = RIIp_Rp8_Br_cAMPS_C_2;
    x_solver[16] = RIIp_Rp8_pCPT_cAMPS_C_2;
    x_solver[17] = RIIp_Rp_cAMPS_C_2;
    x_solver[18] = RIIp_Sp8_Br_cAMPS_C_2;
    x_solver[19] = RIIp_Sp8_Br_cAMPS_2;
    x_solver[20] = RIIp_2;
    x_solver[21] = RII_2;
    x_solver[22] = IBMX;
    x_solver[23] = PDE;
    x_solver[24] = pPDE;
}

} // namespace model_Isensee_JCB2018
} // namespace amici
