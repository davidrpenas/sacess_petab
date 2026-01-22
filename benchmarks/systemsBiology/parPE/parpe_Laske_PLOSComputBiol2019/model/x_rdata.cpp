#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void x_rdata_Laske_PLOSComputBiol2019(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = V_ex;
    x_rdata[1] = V_end;
    x_rdata[2] = Vp_cyt;
    x_rdata[3] = Vp_nuc;
    x_rdata[4] = R_C;
    x_rdata[5] = R_C_RdRp;
    x_rdata[6] = R_V;
    x_rdata[7] = R_V_RdRp;
    x_rdata[8] = R_M1;
    x_rdata[9] = R_M2;
    x_rdata[10] = R_M3;
    x_rdata[11] = R_M4;
    x_rdata[12] = R_M5;
    x_rdata[13] = R_M6;
    x_rdata[14] = R_M7;
    x_rdata[15] = R_M8;
    x_rdata[16] = P_B1;
    x_rdata[17] = P_B2;
    x_rdata[18] = P_PA;
    x_rdata[19] = P_HA;
    x_rdata[20] = P_NP;
    x_rdata[21] = P_NA;
    x_rdata[22] = P_M1;
    x_rdata[23] = P_NEP;
    x_rdata[24] = Vp_nuc_M1;
    x_rdata[25] = Vp_cyt_M1;
    x_rdata[26] = V_rel;
    x_rdata[27] = P_RdRp;
    x_rdata[28] = Cp;
    x_rdata[29] = V_att_Hi;
    x_rdata[30] = V_att_Lo;
    x_rdata[31] = P_M2;
    x_rdata[32] = B_att_Lo;
    x_rdata[33] = B_att_Hi;
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
