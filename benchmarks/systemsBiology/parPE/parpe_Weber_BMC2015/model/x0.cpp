#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Weber_BMC2015 {

void x0_Weber_BMC2015(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 466534.79940000002;
    x0[1] = 123.8608;
    x0[2] = 1577540.5393999999;
    x0[3] = 332054.50410000002;
    x0[4] = 31948388.5902;
    x0[5] = 160797.73639999999;
    x0[6] = 42082828.668099999;
}

} // namespace model_Weber_BMC2015
} // namespace amici
