#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void y_Beer_MolBioSystems2014(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = Bac;
    y[1] = Ind;
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
