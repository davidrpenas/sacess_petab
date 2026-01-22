#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void root_Smith_BMCSystBiol2013(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = t - t_ins;
    root[1] = -PI345P3 + pip3_basal;
    root[2] = PI345P3 - pip3_basal;
    root[3] = t - 2880;
    root[4] = t - 2895;
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
