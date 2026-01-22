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
namespace model_Zhao_QuantBiol2020 {

void y_Zhao_QuantBiol2020(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Confirmed_Infected;
    y[1] = Susceptible;
    y[2] = Unquarantined_Infected;
    y[3] = Quarantined_Infected;
    y[4] = Cumulative_Infected;
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
