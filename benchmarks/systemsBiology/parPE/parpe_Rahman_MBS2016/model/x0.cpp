#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Rahman_MBS2016 {

void x0_Rahman_MBS2016(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 17940000.0;
    x0[1] = 16300.0;
    x0[2] = 9000.0;
    x0[3] = 11000.0;
}

} // namespace model_Rahman_MBS2016
} // namespace amici
