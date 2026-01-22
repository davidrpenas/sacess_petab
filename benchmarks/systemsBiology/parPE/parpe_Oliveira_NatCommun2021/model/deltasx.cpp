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
#include "xdot.h"
#include "xdot_old.h"
#include "sx.h"
#include "stau.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void deltasx_Oliveira_NatCommun2021(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl){
    switch(ie) {
        case 0:
        case 1:
            switch(ip) {
                case 7:
                    deltasx[0] = -dExposeddt + xdot_old0;
                    deltasx[1] = -dSusceptibledt + xdot_old1;
                    deltasx[2] = -dAsymptomaticdt + xdot_old2;
                    deltasx[3] = -dSymptomaticdt + xdot_old3;
                    deltasx[4] = -dRecovereddt + xdot_old4;
                    deltasx[5] = -dHospitaldt + xdot_old5;
                    deltasx[6] = -dICUdt + xdot_old6;
                    deltasx[7] = -dDeathsdt + xdot_old7;
                    deltasx[8] = -dCumulative_casesdt + xdot_old8;
                    break;
            }
            break;
        case 2:
        case 3:
            switch(ip) {
                case 8:
                    deltasx[0] = -dExposeddt + xdot_old0;
                    deltasx[1] = -dSusceptibledt + xdot_old1;
                    deltasx[2] = -dAsymptomaticdt + xdot_old2;
                    deltasx[3] = -dSymptomaticdt + xdot_old3;
                    deltasx[4] = -dRecovereddt + xdot_old4;
                    deltasx[5] = -dHospitaldt + xdot_old5;
                    deltasx[6] = -dICUdt + xdot_old6;
                    deltasx[7] = -dDeathsdt + xdot_old7;
                    deltasx[8] = -dCumulative_casesdt + xdot_old8;
                    break;
            }
            break;
    }
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
