#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "amici/splinefunctions.h"
#include <vector>
#include "p.h"

namespace amici {
namespace model_Crauste_CellSystems2017 {

std::vector<HermiteSpline> create_splines_Crauste_CellSystems2017(const realtype *p, const realtype *k){
    return {};
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
