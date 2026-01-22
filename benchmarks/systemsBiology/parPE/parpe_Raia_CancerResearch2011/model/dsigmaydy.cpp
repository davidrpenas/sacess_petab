#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"

namespace amici {
namespace model_Raia_CancerResearch2011 {

void dsigmaydy_Raia_CancerResearch2011(realtype *dsigmaydy, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    dsigmaydy[0] = noiseParameter2_observable_CD274mRNA;
    dsigmaydy[9] = noiseParameter2_observable_IL13_cell;
    dsigmaydy[18] = noiseParameter2_observable_RecSurf;
    dsigmaydy[27] = noiseParameter2_observable_SOCS3;
    dsigmaydy[36] = noiseParameter2_observable_SOCS3mRNA;
    dsigmaydy[45] = noiseParameter2_observable_pIL4Ra;
    dsigmaydy[54] = noiseParameter2_observable_pJAK2;
    dsigmaydy[63] = noiseParameter2_observable_pSTAT5;
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
