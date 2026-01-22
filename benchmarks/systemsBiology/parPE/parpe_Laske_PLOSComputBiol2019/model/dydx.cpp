#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

void dydx_Laske_PLOSComputBiol2019(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[14] = 1;
    dydx[18] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[19] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[20] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[21] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[22] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[23] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[24] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[25] = -800*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[27] = 1.0/8.0;
    dydx[31] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[32] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[33] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[34] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[35] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[36] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[37] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[38] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[40] = 1.0/8.0;
    dydx[44] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[45] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[46] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[47] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[48] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[49] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[50] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[51] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[54] = 1.0/8.0;
    dydx[67] = 1.0/8.0;
    dydx[79] = 1.0/8.0;
    dydx[156] = 1;
    dydx[313] = 1.0/8.0;
    dydx[317] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[318] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[319] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[320] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[321] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[322] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[323] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[324] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2) + 100/(RNP_cyt + RNP_nuc);
    dydx[326] = 1.0/8.0;
    dydx[330] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[331] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[332] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[333] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[334] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[335] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[336] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[337] = -100*RNP_nuc/std::pow(RNP_cyt + RNP_nuc, 2);
    dydx[341] = 1;
    dydx[366] = 1.0/8.0;
    dydx[378] = 1;
    dydx[391] = 1;
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
