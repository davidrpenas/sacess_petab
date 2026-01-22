#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void y_Bruno_JExpBot2016(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = b10;
    y[1] = b10;
    y[2] = b10;
    y[3] = b10;
    y[4] = b10;
    y[5] = b10;
    y[6] = b10;
    y[7] = b10;
    y[8] = b10;
    y[9] = b10;
    y[10] = b10;
    y[11] = b10;
    y[12] = b10;
    y[13] = b10;
    y[14] = b10;
    y[15] = b10;
    y[16] = b10;
    y[17] = b10;
    y[18] = b10;
    y[19] = b10;
    y[20] = b10;
    y[21] = b10;
    y[22] = b10;
    y[23] = b10;
    y[24] = b10;
    y[25] = b10;
    y[26] = b10;
    y[27] = b10;
    y[28] = bcar;
    y[29] = bcar;
    y[30] = bcar;
    y[31] = bcar;
    y[32] = bcar;
    y[33] = bcar;
    y[34] = bcar;
    y[35] = bcar;
    y[36] = bcar;
    y[37] = bcar;
    y[38] = bcar;
    y[39] = bcar;
    y[40] = bcar;
    y[41] = bcar;
    y[42] = bcry;
    y[43] = bcry;
    y[44] = bcry;
    y[45] = bcry;
    y[46] = bcry;
    y[47] = bcry;
    y[48] = bcry;
    y[49] = ohb10;
    y[50] = ohb10;
    y[51] = ohb10;
    y[52] = ohb10;
    y[53] = ohb10;
    y[54] = ohb10;
    y[55] = ohb10;
    y[56] = ohb10;
    y[57] = ohb10;
    y[58] = ohb10;
    y[59] = ohb10;
    y[60] = ohb10;
    y[61] = ohb10;
    y[62] = ohb10;
    y[63] = ohb10;
    y[64] = ohb10;
    y[65] = ohb10;
    y[66] = ohb10;
    y[67] = ohb10;
    y[68] = ohb10;
    y[69] = ohb10;
    y[70] = zea;
    y[71] = zea;
    y[72] = zea;
    y[73] = zea;
    y[74] = zea;
    y[75] = zea;
    y[76] = zea;
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
