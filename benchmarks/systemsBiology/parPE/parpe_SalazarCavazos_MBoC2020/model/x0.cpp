#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void x0_SalazarCavazos_MBoC2020(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = EGFconc*NA*Vextra;
    x0[1] = EGFR_total;
    x0[2] = GRB2_f*GRB2_total__FREE*f;
    x0[3] = SHC1_total__FREE*f;
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
