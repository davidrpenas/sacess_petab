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
namespace model_Schwen_PONE2014 {

void dydx_Schwen_PONE2014(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[3] = observableParameter2_observable_Insulin/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin) - observableParameter2_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(observableParameter4_observable_Insulin*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
    dydx[30] = observableParameter2_observable_IR1;
    dydx[32] = 0.60499999999999998*observableParameter2_observable_IRsum;
    dydx[41] = observableParameter2_observable_IR2;
    dydx[42] = 0.39500000000000002*observableParameter2_observable_IRsum;
    dydx[50] = observableParameter2_observable_IR1;
    dydx[52] = 0.60499999999999998*observableParameter2_observable_IRsum;
    dydx[61] = observableParameter2_observable_IR2;
    dydx[62] = 0.39500000000000002*observableParameter2_observable_IRsum;
    dydx[93] = observableParameter2_observable_Insulin*observableParameter3_observable_Insulin/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin) - observableParameter2_observable_Insulin*observableParameter3_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(observableParameter4_observable_Insulin*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
}

} // namespace model_Schwen_PONE2014
} // namespace amici
