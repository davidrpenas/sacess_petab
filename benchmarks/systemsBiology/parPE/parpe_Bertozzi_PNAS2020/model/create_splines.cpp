#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "amici/splinefunctions.h"
#include <vector>
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Bertozzi_PNAS2020 {

std::vector<HermiteSpline> create_splines_Bertozzi_PNAS2020(const realtype *p, const realtype *k){
    return {};
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
