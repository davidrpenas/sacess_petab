#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void y_Elowitz_Nature2000(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = GFP*observableParameter2_fluorescence + observableParameter1_fluorescence;
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
