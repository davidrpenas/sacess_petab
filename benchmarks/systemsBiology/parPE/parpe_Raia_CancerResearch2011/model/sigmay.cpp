#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void sigmay_Raia_CancerResearch2011(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_observable_CD274mRNA = noiseParameter1_observable_CD274mRNA + noiseParameter2_observable_CD274mRNA*observable_CD274mRNA;  // sigmay[0]
    sigma_observable_IL13_cell = noiseParameter1_observable_IL13_cell + noiseParameter2_observable_IL13_cell*observable_IL13_cell;  // sigmay[1]
    sigma_observable_RecSurf = noiseParameter1_observable_RecSurf + noiseParameter2_observable_RecSurf*observable_RecSurf;  // sigmay[2]
    sigma_observable_SOCS3 = noiseParameter1_observable_SOCS3 + noiseParameter2_observable_SOCS3*observable_SOCS3;  // sigmay[3]
    sigma_observable_SOCS3mRNA = noiseParameter1_observable_SOCS3mRNA + noiseParameter2_observable_SOCS3mRNA*observable_SOCS3mRNA;  // sigmay[4]
    sigma_observable_pIL4Ra = noiseParameter1_observable_pIL4Ra + noiseParameter2_observable_pIL4Ra*observable_pIL4Ra;  // sigmay[5]
    sigma_observable_pJAK2 = noiseParameter1_observable_pJAK2 + noiseParameter2_observable_pJAK2*observable_pJAK2;  // sigmay[6]
    sigma_observable_pSTAT5 = noiseParameter1_observable_pSTAT5 + noiseParameter2_observable_pSTAT5*observable_pSTAT5;  // sigmay[7]
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
