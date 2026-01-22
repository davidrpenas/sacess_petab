#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Schwen_PONE2014 {

void x0_Schwen_PONE2014(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_Ins;
    x0[1] = ini_R1;
    x0[2] = ini_R1*ini_R2fold;
}

} // namespace model_Schwen_PONE2014
} // namespace amici
