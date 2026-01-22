#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"

namespace amici {
namespace model_Perelson_Science1996 {

void w_Perelson_Science1996(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    flux_v1 = K0*T0*Vin;  // w[0]
    flux_v2 = Tstar*delta;  // w[1]
    flux_v3 = Vin*c;  // w[2]
    flux_v4 = Vni*c;  // w[3]
    flux_v5 = NN*Tstar*delta;  // w[4]
}

} // namespace model_Perelson_Science1996
} // namespace amici
