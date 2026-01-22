#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void x0_Beer_MolBioSystems2014(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = init_Bac;
    x0[1] = 10.0;
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
