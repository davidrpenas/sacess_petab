#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdx.h"

namespace amici {
namespace model_Rahman_MBS2016 {

void dydx_Rahman_MBS2016(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2) - 100/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);
    dydx[1] = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);
    dydx[2] = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);
    dydx[3] = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);
    dydx[4] = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);
    dydx[5] = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);
    dydx[6] = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);
}

} // namespace model_Rahman_MBS2016
} // namespace amici
