#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "xdot.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void deltax_Smith_BMCSystBiol2013(double *deltax, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ie, const realtype *xdot, const realtype *xdot_old){
    switch(ie) {
        case 0:
        case 4:
            deltax[1] = -Ins;
            break;
        case 3:
            deltax[1] = 499999 - Ins;
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
