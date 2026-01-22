#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Lucarelli_CellSystems2018 {

void sx0_Lucarelli_CellSystems2018(realtype *sx0, const realtype t, const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 1:
            sx0[3] = 1;
            break;
        case 2:
            sx0[4] = 1;
            break;
        case 3:
            sx0[5] = 1;
            break;
        case 67:
            sx0[1] = 1;
            break;
    }
}

} // namespace model_Lucarelli_CellSystems2018
} // namespace amici
