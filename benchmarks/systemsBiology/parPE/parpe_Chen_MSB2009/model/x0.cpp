#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Chen_MSB2009 {

void x0_Chen_MSB2009(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[13] = 462000.0;
    x0[15] = 6230.0;
    x0[17] = 794.0;
    x0[34] = 1200000000.0;
    x0[72] = 4498.7299999999996;
    x0[276] = 534751.0;
    x0[289] = 1264.9100000000001;
    x0[299] = 58095.199999999997;
    x0[304] = 1100000.0;
    x0[311] = 71131.199999999997;
    x0[317] = 88914000.0;
    x0[318] = 50000.0;
    x0[323] = 3020000.0;
    x0[330] = 124480.0;
    x0[335] = 695000.0;
    x0[344] = 16870.200000000001;
    x0[367] = 35565600.0;
    x0[394] = 905000.0;
    x0[398] = 300416000.0;
    x0[403] = 450000.0;
    x0[406] = 70000.0;
    x0[415] = 1080000.0;
    x0[423] = 2213.5900000000001;
    x0[424] = 393639.0;
    x0[426] = 56100.900000000001;
    x0[428] = 94868.300000000003;
    x0[430] = 1000000.0;
}

} // namespace model_Chen_MSB2009
} // namespace amici
