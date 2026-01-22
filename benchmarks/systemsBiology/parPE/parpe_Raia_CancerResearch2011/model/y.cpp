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

void y_Raia_CancerResearch2011(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = CD274mRNA*observableParameter1_observable_CD274mRNA;
    y[1] = observableParameter1_observable_IL13_cell*(IL13_DecoyR + IL13_Rec + p_IL13_Rec + p_IL13_Rec_i);
    y[2] = IL13_Rec + Rec + p_IL13_Rec;
    y[3] = SOCS3;
    y[4] = SOCS3mRNA*observableParameter1_observable_SOCS3mRNA;
    y[5] = observableParameter1_observable_pIL4Ra*(p_IL13_Rec + p_IL13_Rec_i);
    y[6] = observableParameter1_observable_pJAK2*pJAK2;
    y[7] = pSTAT5;
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
