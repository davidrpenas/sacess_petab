#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void sx0_fixedParameters_SalazarCavazos_MBoC2020(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs){
    static const std::array<int, 4> _x0_fixedParameters_idxs = {
        0, 1, 2, 3
    };
    for(auto idx: reinitialization_state_idxs) {
        if(std::find(_x0_fixedParameters_idxs.cbegin(), _x0_fixedParameters_idxs.cend(), idx) != _x0_fixedParameters_idxs.cend())
            sx0_fixedParameters[idx] = 0.0;
    }
    switch(ip) {
        case 0:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[2] = GRB2_f*f;
            break;
        case 1:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 3) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[3] = f;
            break;
    }
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
