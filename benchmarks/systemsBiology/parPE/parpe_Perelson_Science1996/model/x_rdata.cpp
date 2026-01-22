#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Perelson_Science1996 {

void x_rdata_Perelson_Science1996(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Tstar;
    x_rdata[1] = V;
    x_rdata[2] = Vin;
    x_rdata[3] = Vni;
}

} // namespace model_Perelson_Science1996
} // namespace amici
