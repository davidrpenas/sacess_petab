#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void x_solver_Laske_PLOSComputBiol2019(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = V_ex;
    x_solver[1] = V_end;
    x_solver[2] = Vp_cyt;
    x_solver[3] = Vp_nuc;
    x_solver[4] = R_C;
    x_solver[5] = R_C_RdRp;
    x_solver[6] = R_V;
    x_solver[7] = R_V_RdRp;
    x_solver[8] = R_M1;
    x_solver[9] = R_M2;
    x_solver[10] = R_M3;
    x_solver[11] = R_M4;
    x_solver[12] = R_M5;
    x_solver[13] = R_M6;
    x_solver[14] = R_M7;
    x_solver[15] = R_M8;
    x_solver[16] = P_B1;
    x_solver[17] = P_B2;
    x_solver[18] = P_PA;
    x_solver[19] = P_HA;
    x_solver[20] = P_NP;
    x_solver[21] = P_NA;
    x_solver[22] = P_M1;
    x_solver[23] = P_NEP;
    x_solver[24] = Vp_nuc_M1;
    x_solver[25] = Vp_cyt_M1;
    x_solver[26] = V_rel;
    x_solver[27] = P_RdRp;
    x_solver[28] = Cp;
    x_solver[29] = V_att_Hi;
    x_solver[30] = V_att_Lo;
    x_solver[31] = P_M2;
    x_solver[32] = B_att_Lo;
    x_solver[33] = B_att_Hi;
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
