#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void x0_Oliveira_NatCommun2021(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = exposed_init_concentration*population;
    x0[1] = population*(-asymptomatic_init_concentration - exposed_init_concentration - symptomatic_init_concentration + 1);
    x0[2] = asymptomatic_init_concentration*population;
    x0[3] = population*symptomatic_init_concentration;
    x0[8] = population*symptomatic_init_concentration;
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
