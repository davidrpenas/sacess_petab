#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void x_solver_Lucarelli_CellSystems2018(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = TGFb;
    x_solver[1] = Rec;
    x_solver[2] = TGFb_pRec;
    x_solver[3] = S2;
    x_solver[4] = S3;
    x_solver[5] = S4;
    x_solver[6] = S2_S4_S4;
    x_solver[7] = ppS2_ppS2_ppS2;
    x_solver[8] = ppS3_ppS3_ppS3;
    x_solver[9] = S4_S4_S4;
    x_solver[10] = pS2;
    x_solver[11] = pS3;
    x_solver[12] = ppS2;
    x_solver[13] = ppS3;
    x_solver[14] = ppS2_ppS2_S4;
    x_solver[15] = ppS2_ppS2_ppS3;
    x_solver[16] = ppS2_ppS3_ppS3;
    x_solver[17] = ppS3_ppS3_S4;
    x_solver[18] = ppS2_ppS3_S4;
    x_solver[19] = ppS3_S4_S4;
    x_solver[20] = ppS2_S4_S4;
    x_solver[21] = geneA;
    x_solver[22] = geneB;
    x_solver[23] = geneC;
    x_solver[24] = geneD;
    x_solver[25] = geneE;
    x_solver[26] = geneF;
    x_solver[27] = geneG;
    x_solver[28] = geneH;
    x_solver[29] = geneI;
    x_solver[30] = geneJ;
    x_solver[31] = geneK;
    x_solver[32] = geneL;
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
