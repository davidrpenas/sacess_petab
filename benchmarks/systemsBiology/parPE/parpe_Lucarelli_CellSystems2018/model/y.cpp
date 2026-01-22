#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void y_Lucarelli_CellSystems2018(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = geneH + 9.9999999999999995e-7;
    y[1] = geneI + 9.9999999999999995e-7;
    y[2] = geneC + 9.9999999999999995e-7;
    y[3] = geneJ + 9.9999999999999995e-7;
    y[4] = geneE + 9.9999999999999995e-7;
    y[5] = geneG + 9.9999999999999995e-7;
    y[6] = geneL + 9.9999999999999995e-7;
    y[7] = (100*S2 + 100*S2_S4_S4)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[8] = (100*S2 + 100*S2_S4_S4)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[9] = 100*S3/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[10] = 100*S3/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[11] = 100*pS2/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[12] = 100*pS2/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[13] = 100*pS3/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[14] = 100*pS3/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[15] = (100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[16] = (100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[17] = (100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3 + 100*ppS3 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4 + 300*ppS3_ppS3_ppS3)/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[18] = (100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3 + 100*ppS3 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4 + 300*ppS3_ppS3_ppS3)/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[19] = S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3;
    y[20] = (100*S2 + 100*S2_S4_S4 + 100*pS2 + 100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 2*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[21] = (100*S2 + 100*S2_S4_S4 + 100*pS2 + 100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + 3*S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 4*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[22] = (100*S3 + 100*pS3 + 100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3 + 100*ppS3 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4 + 300*ppS3_ppS3_ppS3)/(S2 + S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 2*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[23] = (100*S3 + 100*pS3 + 100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3 + 100*ppS3 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4 + 300*ppS3_ppS3_ppS3)/(S2 + 3*S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 4*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[24] = (200*S2_S4_S4 + 200*ppS2_S4_S4 + 100*ppS2_ppS2_S4 + 200*ppS2_ppS3_S4 + 200*ppS3_S4_S4 + 100*ppS3_ppS3_S4)/(S2 + 3*S2_S4_S4 + S3 + pS2 + pS3 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 4*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[25] = (100*S2 + 100*S2_S4_S4)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[26] = (100*S2 + 100*S2_S4_S4)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[28] = 100*pS2/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[29] = 100*pS2/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[31] = (100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[32] = (100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + S2_S4_S4 + pS2 + ppS2 + ppS2_S4_S4 + 2*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 2*ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + ppS2_ppS3_ppS3);
    y[33] = 100.0;
    y[34] = (100*S2 + 100*S2_S4_S4 + 100*pS2 + 100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + 3*S2_S4_S4 + pS2 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3);
    y[35] = (100*S2 + 100*S2_S4_S4 + 100*pS2 + 100*ppS2 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 300*ppS2_ppS2_ppS2 + 200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S2 + 3*S2_S4_S4 + pS2 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3);
    y[36] = (100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3)/(S2 + 3*S2_S4_S4 + pS2 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3);
    y[37] = (100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3)/(S2 + 3*S2_S4_S4 + pS2 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3);
    y[38] = (200*ppS2_S4_S4 + 100*ppS2_ppS2_S4 + 100*ppS2_ppS3_S4)/(S2 + 3*S2_S4_S4 + pS2 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3);
    y[39] = (200*ppS2_S4_S4 + 100*ppS2_ppS2_S4 + 100*ppS2_ppS3_S4)/(S2 + 3*S2_S4_S4 + pS2 + ppS2 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS2_ppS2 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3);
    y[41] = 100*S3/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[43] = 100*pS3/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[44] = 100.0;
    y[45] = (100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3 + 100*ppS3 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4 + 300*ppS3_ppS3_ppS3)/(S3 + pS3 + ppS2_ppS2_ppS3 + ppS2_ppS3_S4 + 2*ppS2_ppS3_ppS3 + ppS3 + ppS3_S4_S4 + 2*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[46] = (200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S3 + pS3 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[47] = (200*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 100*ppS2_ppS3_ppS3)/(S3 + pS3 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[48] = (100*S3 + 100*pS3 + 100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3 + 100*ppS3 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4 + 300*ppS3_ppS3_ppS3)/(S3 + pS3 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[49] = (100*S3 + 100*pS3 + 100*ppS2_ppS2_ppS3 + 100*ppS2_ppS3_S4 + 200*ppS2_ppS3_ppS3 + 100*ppS3 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4 + 300*ppS3_ppS3_ppS3)/(S3 + pS3 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[50] = (100*ppS2_ppS3_S4 + 200*ppS3_S4_S4 + 100*ppS3_ppS3_S4)/(S3 + pS3 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[51] = (100*ppS2_ppS3_S4 + 200*ppS3_S4_S4 + 100*ppS3_ppS3_S4)/(S3 + pS3 + 3*ppS2_ppS2_ppS3 + 3*ppS2_ppS3_S4 + 3*ppS2_ppS3_ppS3 + ppS3 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4 + 3*ppS3_ppS3_ppS3);
    y[52] = (100*S2_S4_S4 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 100*ppS2_ppS3_S4)/(3*S2_S4_S4 + S4 + 3*S4_S4_S4 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS3_S4 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4);
    y[53] = (100*S2_S4_S4 + 100*ppS2_S4_S4 + 200*ppS2_ppS2_S4 + 100*ppS2_ppS3_S4)/(3*S2_S4_S4 + S4 + 3*S4_S4_S4 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS3_S4 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4);
    y[54] = (100*ppS2_ppS3_S4 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4)/(3*S2_S4_S4 + S4 + 3*S4_S4_S4 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS3_S4 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4);
    y[55] = (100*ppS2_ppS3_S4 + 100*ppS3_S4_S4 + 200*ppS3_ppS3_S4)/(3*S2_S4_S4 + S4 + 3*S4_S4_S4 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS3_S4 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4);
    y[56] = 2*S2_S4_S4 + S4 + 3*S4_S4_S4 + 2*ppS2_S4_S4 + ppS2_ppS2_S4 + ppS2_ppS3_S4 + 2*ppS3_S4_S4 + ppS3_ppS3_S4;
    y[57] = (200*S2_S4_S4 + 100*S4 + 300*S4_S4_S4 + 200*ppS2_S4_S4 + 100*ppS2_ppS2_S4 + 100*ppS2_ppS3_S4 + 200*ppS3_S4_S4 + 100*ppS3_ppS3_S4)/(3*S2_S4_S4 + S4 + 3*S4_S4_S4 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS3_S4 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4);
    y[58] = (200*S2_S4_S4 + 100*S4 + 300*S4_S4_S4 + 200*ppS2_S4_S4 + 100*ppS2_ppS2_S4 + 100*ppS2_ppS3_S4 + 200*ppS3_S4_S4 + 100*ppS3_ppS3_S4)/(3*S2_S4_S4 + S4 + 3*S4_S4_S4 + 3*ppS2_S4_S4 + 3*ppS2_ppS2_S4 + 3*ppS2_ppS3_S4 + 3*ppS3_S4_S4 + 3*ppS3_ppS3_S4);
    y[59] = geneA + 9.9999999999999995e-7;
    y[60] = geneB + 9.9999999999999995e-7;
    y[61] = geneF + 9.9999999999999995e-7;
    y[62] = geneD + 9.9999999999999995e-7;
    y[63] = Rec + TGFb_pRec;
    y[64] = geneK + 9.9999999999999995e-7;
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
