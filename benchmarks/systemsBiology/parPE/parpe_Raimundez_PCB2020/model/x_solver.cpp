#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void x_solver_Raimundez_PCB2020(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EGFR;
    x_solver[1] = EGFR_CET;
    x_solver[2] = EGFR_EGF;
    x_solver[3] = EGFR_EGF_2;
    x_solver[4] = pEGFR_EGF_2;
    x_solver[5] = pEGFR_EGF_2_i;
    x_solver[6] = MMET;
    x_solver[7] = MMET_2;
    x_solver[8] = pMMET_2;
    x_solver[9] = pMMET_2_i;
    x_solver[10] = MMET_EGFR;
    x_solver[11] = pMMET_pEGFR;
    x_solver[12] = pMMET_pEGFR_i;
    x_solver[13] = MMET_METinh;
    x_solver[14] = MMET_MMET_METinh;
    x_solver[15] = MMET_METinh_2;
    x_solver[16] = EGFR_MMET_METinh;
    x_solver[17] = RAS_GTP;
    x_solver[18] = pMAPK;
    x_solver[19] = pPI3K;
    x_solver[20] = pMPI3K;
    x_solver[21] = pAKT;
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
