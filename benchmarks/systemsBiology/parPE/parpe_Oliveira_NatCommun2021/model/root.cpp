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
namespace model_Oliveira_NatCommun2021 {

void root_Oliveira_NatCommun2021(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    root[0] = t - t_1;
    root[1] = -t + t_1;
    root[2] = t - t_2;
    root[3] = -t + t_2;
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
