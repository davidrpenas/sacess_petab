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
namespace model_Sneyd_PNAS2002 {

void y_Sneyd_PNAS2002(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = 0.65610000000000002*std::pow(IPR_A + 0.11111111111111112*IPR_O, 4);
}

} // namespace model_Sneyd_PNAS2002
} // namespace amici
