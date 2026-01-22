#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void y_Borghans_BiophysChem1997(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Z_state*observableParameter2_Ca + observableParameter1_Ca;
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
