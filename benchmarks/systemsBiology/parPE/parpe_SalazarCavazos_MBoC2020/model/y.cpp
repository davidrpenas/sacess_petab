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
namespace model_SalazarCavazos_MBoC2020 {

void y_SalazarCavazos_MBoC2020(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = 100*pY1068/EGFRtot;
    y[1] = 100*pY1068/EGFRtot;
    y[2] = 100*pY1068/EGFRtot;
    y[3] = 100*pY1068/EGFRtot;
    y[4] = 100*pY1068/EGFRtot;
    y[5] = 100*pY1068/EGFRtot;
    y[6] = 100*pY1173/EGFRtot;
    y[7] = 100*pY1173/EGFRtot;
    y[8] = 100*pY1173/EGFRtot;
    y[9] = 100*pY1173/EGFRtot;
    y[10] = 100*pY1173/EGFRtot;
    y[11] = 100*pY1173/EGFRtot;
    y[12] = (100*EGFRtot - 100*unphosR)/EGFRtot;
    y[13] = (100*EGFRtot - 100*unphosR)/EGFRtot;
    y[14] = (100*EGFRtot - 100*unphosR)/EGFRtot;
    y[15] = (100*EGFRtot - 100*unphosR)/EGFRtot;
    y[16] = (100*EGFRtot - 100*unphosR)/EGFRtot;
    y[17] = (100*EGFRtot - 100*unphosR)/EGFRtot;
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
