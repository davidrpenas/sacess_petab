#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void x_rdata_Lucarelli_CellSystems2018(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = TGFb;
    x_rdata[1] = Rec;
    x_rdata[2] = TGFb_pRec;
    x_rdata[3] = S2;
    x_rdata[4] = S3;
    x_rdata[5] = S4;
    x_rdata[6] = S2_S4_S4;
    x_rdata[7] = ppS2_ppS2_ppS2;
    x_rdata[8] = ppS3_ppS3_ppS3;
    x_rdata[9] = S4_S4_S4;
    x_rdata[10] = pS2;
    x_rdata[11] = pS3;
    x_rdata[12] = ppS2;
    x_rdata[13] = ppS3;
    x_rdata[14] = ppS2_ppS2_S4;
    x_rdata[15] = ppS2_ppS2_ppS3;
    x_rdata[16] = ppS2_ppS3_ppS3;
    x_rdata[17] = ppS3_ppS3_S4;
    x_rdata[18] = ppS2_ppS3_S4;
    x_rdata[19] = ppS3_S4_S4;
    x_rdata[20] = ppS2_S4_S4;
    x_rdata[21] = geneA;
    x_rdata[22] = geneB;
    x_rdata[23] = geneC;
    x_rdata[24] = geneD;
    x_rdata[25] = geneE;
    x_rdata[26] = geneF;
    x_rdata[27] = geneG;
    x_rdata[28] = geneH;
    x_rdata[29] = geneI;
    x_rdata[30] = geneJ;
    x_rdata[31] = geneK;
    x_rdata[32] = geneL;
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
