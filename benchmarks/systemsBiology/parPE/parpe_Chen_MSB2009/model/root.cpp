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
namespace model_Chen_MSB2009 {

void root_Chen_MSB2009(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = t - 1800;
    root[1] = 1800 - t;
    root[2] = t - 1799.99;
    root[3] = 1799.99 - t;
    root[4] = t - 2660;
    root[5] = 2660 - t;
    root[6] = t - 2659.9899999999998;
    root[7] = 2659.9899999999998 - t;
}

} // namespace model_Chen_MSB2009
} // namespace amici
