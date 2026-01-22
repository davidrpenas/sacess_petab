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
namespace model_Weber_BMC2015 {

void root_Weber_BMC2015(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = -PdBu_time + t;
    root[1] = PdBu_time - t;
    root[2] = -kb_NB142_70_time + t;
    root[3] = kb_NB142_70_time - t;
    root[4] = t;
    root[5] = -t;
}

} // namespace model_Weber_BMC2015
} // namespace amici
