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

namespace amici {
namespace model_Weber_BMC2015 {

void y_Weber_BMC2015(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CERT*observableParameter1_yCERTpRN24/(CERT + CERTERa + CERTTGNa);
    y[1] = CERT + CERTERa + CERTTGNa;
    y[2] = PI4K3Ba*observableParameter1_yPI4K3BpRN24/(PI4K3B + PI4K3Ba);
    y[3] = PI4K3B + PI4K3Ba;
    y[4] = PKDDAGa*observableParameter1_yPKDpN0;
    y[5] = PKDDAGa*observableParameter1_yPKDpN24;
    y[6] = PKDDAGa*observableParameter1_yPKDpN25;
    y[7] = PKD + PKDDAGa;
}

} // namespace model_Weber_BMC2015
} // namespace amici
