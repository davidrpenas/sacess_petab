#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Raimundez_PCB2020 {

void x_rdata_Raimundez_PCB2020(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = EGFR;
    x_rdata[1] = EGFR_CET;
    x_rdata[2] = EGFR_EGF;
    x_rdata[3] = EGFR_EGF_2;
    x_rdata[4] = pEGFR_EGF_2;
    x_rdata[5] = pEGFR_EGF_2_i;
    x_rdata[6] = MMET;
    x_rdata[7] = MMET_2;
    x_rdata[8] = pMMET_2;
    x_rdata[9] = pMMET_2_i;
    x_rdata[10] = MMET_EGFR;
    x_rdata[11] = pMMET_pEGFR;
    x_rdata[12] = pMMET_pEGFR_i;
    x_rdata[13] = MMET_METinh;
    x_rdata[14] = MMET_MMET_METinh;
    x_rdata[15] = MMET_METinh_2;
    x_rdata[16] = EGFR_MMET_METinh;
    x_rdata[17] = RAS_GTP;
    x_rdata[18] = pMAPK;
    x_rdata[19] = pPI3K;
    x_rdata[20] = pMPI3K;
    x_rdata[21] = pAKT;
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
