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
namespace model_Giordano_Nature2020 {

void y_Giordano_Nature2020(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CumulativeDiagnosed;
    y[1] = Extinct;
    y[2] = DiagnosedHealed;
    y[3] = CurrentDiagnosedInfected;
    y[4] = Diagnosed;
    y[5] = Recognized;
    y[6] = Threatened;
}

} // namespace model_Giordano_Nature2020
} // namespace amici
