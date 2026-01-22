#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "w.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void dydp_Borghans_BiophysChem1997(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 20:
            dydp[0] = 1;
            break;
        case 21:
            dydp[0] = Z_state;
            break;
    }
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
