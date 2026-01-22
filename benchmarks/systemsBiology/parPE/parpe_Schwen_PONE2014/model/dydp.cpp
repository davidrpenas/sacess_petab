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

void dydp_Schwen_PONE2014(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 0:
            dydp[4] = 1/(ini_R1*M_LN10);
            break;
        case 1:
            dydp[5] = 1/(ini_R2fold*M_LN10);
            break;
        case 2:
            dydp[6] = 1/(ka1*M_LN10);
            break;
        case 3:
            dydp[7] = 1/(ka2fold*M_LN10);
            break;
        case 4:
            dydp[8] = 1/(kd1*M_LN10);
            break;
        case 5:
            dydp[9] = 1/(kd2fold*M_LN10);
            break;
        case 13:
            dydp[0] = observableParameter2_observable_IR1;
            break;
        case 14:
            dydp[0] = IR1 + IR1in + observableParameter1_observable_IR1;
            break;
        case 15:
            dydp[1] = observableParameter2_observable_IR2;
            break;
        case 16:
            dydp[1] = IR2 + IR2in + observableParameter1_observable_IR2;
            break;
        case 17:
            dydp[2] = observableParameter2_observable_IRsum;
            break;
        case 18:
            dydp[2] = 0.60499999999999998*IR1 + 0.60499999999999998*IR1in + 0.39500000000000002*IR2 + 0.39500000000000002*IR2in + observableParameter1_observable_IRsum;
            break;
        case 19:
            dydp[3] = 1;
            break;
        case 20:
            dydp[3] = (Ins + InsulinFragments*observableParameter3_observable_Insulin)/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin);
            break;
        case 21:
            dydp[3] = InsulinFragments*observableParameter2_observable_Insulin/(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin) - InsulinFragments*observableParameter2_observable_Insulin*(Ins + InsulinFragments*observableParameter3_observable_Insulin)/(observableParameter4_observable_Insulin*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
            break;
        case 22:
            dydp[3] = observableParameter2_observable_Insulin*std::pow(Ins + InsulinFragments*observableParameter3_observable_Insulin, 2)/(std::pow(observableParameter4_observable_Insulin, 2)*std::pow(1 + (Ins + InsulinFragments*observableParameter3_observable_Insulin)/observableParameter4_observable_Insulin, 2));
            break;
    }
}

} // namespace model_Schwen_PONE2014
} // namespace amici
