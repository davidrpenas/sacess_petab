#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "w.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void w_Beer_MolBioSystems2014(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){

    // dynamic expressions
    lag = Heaviside_0;  // w[0]
    flux_reaction1 = 1.0*Bac*beta*lag*(-Bac + Bacmax)/Bacmax;  // w[1]
    flux_reaction2 = -1.0*Bac*Glu*ksyn;  // w[2]
    flux_reaction3 = 1.0*Bac*Glu*ksyn - 1.0*std::pow(cGlu, 2)*kdim;  // w[3]
    flux_reaction4 = -1.0*Ind*kdegi + 1.0*std::pow(cGlu, 2)*kdim;  // w[4]
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
