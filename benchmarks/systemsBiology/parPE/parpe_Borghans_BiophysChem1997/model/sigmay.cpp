#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Borghans_BiophysChem1997 {

void sigmay_Borghans_BiophysChem1997(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_Ca = noiseParameter1_Ca;  // sigmay[0]
}

} // namespace model_Borghans_BiophysChem1997
} // namespace amici
