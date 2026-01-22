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
namespace model_Schwen_PONE2014 {

void y_Schwen_PONE2014(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = observableParameter2_observable_IR1*(IR1 + IR1in + observableParameter1_observable_IR1);
    y[1] = observableParameter2_observable_IR2*(IR2 + IR2in + observableParameter1_observable_IR2);
    y[2] = observableParameter2_observable_IRsum*(0.60499999999999998*IR1 + 0.60499999999999998*IR1in + 0.39500000000000002*IR2 + 0.39500000000000002*IR2in + observableParameter1_observable_IRsum);
    y[3] = observableParameter1_observable_Insulin + observableParameter2_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin);
    y[4] = std::log(ini_R1)/M_LN10;
    y[5] = std::log(ini_R2fold)/M_LN10;
    y[6] = std::log(ka1)/M_LN10;
    y[7] = std::log(ka2fold)/M_LN10;
    y[8] = std::log(kd1)/M_LN10;
    y[9] = std::log(kd2fold)/M_LN10;
}

} // namespace model_Schwen_PONE2014
} // namespace amici
