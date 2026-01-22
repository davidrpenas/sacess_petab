#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"

namespace amici {
namespace model_Oliveira_NatCommun2021 {

void sx0_fixedParameters_Oliveira_NatCommun2021(realtype *sx0_fixedParameters, const realtype t, const realtype *x0, const realtype *p, const realtype *k, const int ip, gsl::span<const int> reinitialization_state_idxs){
    static const std::array<int, 5> _x0_fixedParameters_idxs = {
        0, 1, 2, 3, 8
    };
    for(auto idx: reinitialization_state_idxs) {
        if(std::find(_x0_fixedParameters_idxs.cbegin(), _x0_fixedParameters_idxs.cend(), idx) != _x0_fixedParameters_idxs.cend())
            sx0_fixedParameters[idx] = 0.0;
    }
    switch(ip) {
        case 9:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[0] = population;
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[1] = -population;
            break;
        case 10:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[1] = -population;
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[2] = population;
            break;
        case 11:
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[1] = -population;
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 3) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[3] = population;
            if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 8) != reinitialization_state_idxs.cend())
                sx0_fixedParameters[8] = population;
            break;
    }
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
