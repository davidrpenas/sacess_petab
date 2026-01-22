#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Weber_BMC2015 {

void dydx_Weber_BMC2015(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[7] = 1;
    dydx[12] = observableParameter1_yPKDpN0;
    dydx[13] = observableParameter1_yPKDpN24;
    dydx[14] = observableParameter1_yPKDpN25;
    dydx[15] = 1;
    dydx[18] = -PI4K3Ba*observableParameter1_yPI4K3BpRN24/std::pow(PI4K3B + PI4K3Ba, 2);
    dydx[19] = 1;
    dydx[26] = -PI4K3Ba*observableParameter1_yPI4K3BpRN24/std::pow(PI4K3B + PI4K3Ba, 2) + observableParameter1_yPI4K3BpRN24/(PI4K3B + PI4K3Ba);
    dydx[27] = 1;
    dydx[32] = -CERT*observableParameter1_yCERTpRN24/std::pow(CERT + CERTERa + CERTTGNa, 2);
    dydx[33] = 1;
    dydx[40] = -CERT*observableParameter1_yCERTpRN24/std::pow(CERT + CERTERa + CERTTGNa, 2) + observableParameter1_yCERTpRN24/(CERT + CERTERa + CERTTGNa);
    dydx[41] = 1;
    dydx[48] = -CERT*observableParameter1_yCERTpRN24/std::pow(CERT + CERTERa + CERTTGNa, 2);
    dydx[49] = 1;
}

} // namespace model_Weber_BMC2015
} // namespace amici
