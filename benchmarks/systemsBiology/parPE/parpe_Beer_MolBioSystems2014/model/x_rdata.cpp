#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void x_rdata_Beer_MolBioSystems2014(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = Bac;
    x_rdata[1] = Glu;
    x_rdata[2] = cGlu;
    x_rdata[3] = Ind;
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
