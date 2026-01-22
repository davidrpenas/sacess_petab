#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Zheng_PNAS2012 {

void x0_Zheng_PNAS2012(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 0.0041772497634575897;
    x0[1] = 0.010210466858764099;
    x0[2] = 0.0016541281027940701;
    x0[3] = 0.169690316239546;
    x0[4] = 0.0078328187288069003;
    x0[5] = 0.00632744816295157;
    x0[6] = 0.050493521480754401;
    x0[7] = 0.594249755169037;
    x0[8] = 0.0143896310177379;
    x0[9] = 0.0033016841160416499;
    x0[10] = 0.102748675077958;
    x0[11] = 0.026337263499652899;
    x0[12] = 0.0025083103492027701;
    x0[13] = 0.0047183143600213404;
    x0[14] = 0.0013604163179556201;
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
