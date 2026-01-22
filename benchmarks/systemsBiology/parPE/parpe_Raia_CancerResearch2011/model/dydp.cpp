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
namespace model_Raia_CancerResearch2011 {

void dydp_Raia_CancerResearch2011(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 18:
            dydp[0] = CD274mRNA;
            break;
        case 19:
            dydp[1] = IL13_DecoyR + IL13_Rec + p_IL13_Rec + p_IL13_Rec_i;
            break;
        case 20:
            dydp[4] = SOCS3mRNA;
            break;
        case 21:
            dydp[5] = p_IL13_Rec + p_IL13_Rec_i;
            break;
        case 22:
            dydp[6] = pJAK2;
            break;
    }
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
