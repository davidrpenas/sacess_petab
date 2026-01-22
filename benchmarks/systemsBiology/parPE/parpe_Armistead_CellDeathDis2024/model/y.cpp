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
namespace model_Armistead_CellDeathDis2024 {

void y_Armistead_CellDeathDis2024(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Sphinga;
    y[1] = Cer;
    y[2] = Sphingo;
    y[3] = S1P;
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
