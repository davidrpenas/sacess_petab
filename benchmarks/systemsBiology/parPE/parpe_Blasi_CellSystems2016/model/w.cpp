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
namespace model_Blasi_CellSystems2016 {

void w_Blasi_CellSystems2016(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        a_k12 = a_basal;  // w[0]
        a_k12_k12k16 = a_basal;  // w[1]
        a_k12_k8k12 = a_basal;  // w[2]
        a_k12k16_k5k12k16 = a_basal;  // w[3]
        a_k16 = a_basal;  // w[4]
        a_k16_k5k16 = a_basal;  // w[5]
        a_k16_k8k16 = a_basal;  // w[6]
        a_k5 = a_basal;  // w[7]
        a_k5_k5k16 = a_basal;  // w[8]
        a_k5_k5k8 = a_basal;  // w[9]
        a_k5k12_k5k12k16 = a_basal;  // w[10]
        a_k5k12k16_4ac = a_basal;  // w[11]
        a_k5k16_k5k12k16 = a_basal;  // w[12]
        a_k5k16_k5k8k16 = a_basal;  // w[13]
        a_k5k8_k5k8k12 = a_basal;  // w[14]
        a_k5k8_k5k8k16 = a_basal;  // w[15]
        a_k5k8k12_4ac = a_basal;  // w[16]
        a_k5k8k16_4ac = a_basal;  // w[17]
        a_k8_k5k8 = a_basal;  // w[18]
        a_k8_k8k12 = a_basal;  // w[19]
        a_k8_k8k16 = a_basal;  // w[20]
        a_k8k12_k5k8k12 = a_basal;  // w[21]
        a_k8k12_k8k12k16 = a_basal;  // w[22]
        a_k8k16_k5k8k16 = a_basal;  // w[23]
        a_k8k16_k8k12k16 = a_basal;  // w[24]
    }

    // dynamic expressions
    flux_re1 = a_k5*x_0ac - d*x_k5;  // w[25]
    flux_re2 = a_k8*x_0ac - d*x_k8;  // w[26]
    flux_re3 = a_k12*x_0ac - d*x_k12;  // w[27]
    flux_re4 = a_k16*x_0ac - d*x_k16;  // w[28]
    flux_re5 = a_k5_k5k8*x_k5 - d*x_k5k8;  // w[29]
    flux_re6 = a_k5_k5k12*x_k5 - d*x_k5k12;  // w[30]
    flux_re7 = a_k5_k5k16*x_k5 - d*x_k5k16;  // w[31]
    flux_re8 = a_k8_k5k8*x_k8 - d*x_k5k8;  // w[32]
    flux_re9 = a_k8_k8k12*x_k8 - d*x_k8k12;  // w[33]
    flux_re10 = a_k8_k8k16*x_k8 - d*x_k8k16;  // w[34]
    flux_re11 = a_k12_k5k12*x_k12 - d*x_k5k12;  // w[35]
    flux_re12 = a_k12_k8k12*x_k12 - d*x_k8k12;  // w[36]
    flux_re13 = a_k12_k12k16*x_k12 - d*x_k12k16;  // w[37]
    flux_re14 = a_k16_k5k16*x_k16 - d*x_k5k16;  // w[38]
    flux_re15 = a_k16_k8k16*x_k16 - d*x_k8k16;  // w[39]
    flux_re16 = a_k16_k12k16*x_k16 - d*x_k12k16;  // w[40]
    flux_re17 = a_k5k8_k5k8k12*x_k5k8 - d*x_k5k8k12;  // w[41]
    flux_re18 = a_k5k8_k5k8k16*x_k5k8 - d*x_k5k8k16;  // w[42]
    flux_re19 = a_k5k12_k5k8k12*x_k5k12 - d*x_k5k8k12;  // w[43]
    flux_re20 = a_k5k12_k5k12k16*x_k5k12 - d*x_k5k12k16;  // w[44]
    flux_re21 = a_k5k16_k5k8k16*x_k5k16 - d*x_k5k8k16;  // w[45]
    flux_re22 = a_k5k16_k5k12k16*x_k5k16 - d*x_k5k12k16;  // w[46]
    flux_re23 = a_k8k12_k5k8k12*x_k8k12 - d*x_k5k8k12;  // w[47]
    flux_re24 = a_k8k12_k8k12k16*x_k8k12 - d*x_k8k12k16;  // w[48]
    flux_re25 = a_k8k16_k5k8k16*x_k8k16 - d*x_k5k8k16;  // w[49]
    flux_re26 = a_k8k16_k8k12k16*x_k8k16 - d*x_k8k12k16;  // w[50]
    flux_re27 = a_k12k16_k5k12k16*x_k12k16 - d*x_k5k12k16;  // w[51]
    flux_re28 = a_k12k16_k8k12k16*x_k12k16 - d*x_k8k12k16;  // w[52]
    flux_re29 = a_k5k8k12_4ac*x_k5k8k12 - d*x_4ac;  // w[53]
    flux_re30 = a_k5k8k16_4ac*x_k5k8k16 - d*x_4ac;  // w[54]
    flux_re31 = a_k5k12k16_4ac*x_k5k12k16 - d*x_4ac;  // w[55]
    flux_re32 = a_k8k12k16_4ac*x_k8k12k16 - d*x_4ac;  // w[56]
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
