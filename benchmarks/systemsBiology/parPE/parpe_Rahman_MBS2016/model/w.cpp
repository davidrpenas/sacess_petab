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
namespace model_Rahman_MBS2016 {

void w_Rahman_MBS2016(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        infected_normal_transmission_rate = infected_moderate_transmission_rate*infected_normal_transmission_rate_relative;  // w[0]
        infected_weak_transmission_rate = infected_moderate_transmission_rate*infected_weak_transmission_rate_relative;  // w[1]
        treated_transmission_rate = 0.040000000000000001*infected_moderate_transmission_rate;  // w[3]
        flux_v1 = recruitment_rate;  // w[5]
    }

    // dynamic expressions
    prevalence = -100*susceptible/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + 100;  // w[2]
    force_of_infection = (infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);  // w[4]
    flux_v2 = force_of_infection*susceptible;  // w[6]
    flux_v3 = infected_normal*infected_normal_worsen_rate;  // w[7]
    flux_v4 = infected_moderate*infected_moderate_worsen_rate;  // w[8]
    flux_v5 = treated_weak*treated_weak_improve_rate;  // w[9]
    flux_v6 = treated_moderate*treated_moderate_improve_rate;  // w[10]
    flux_v7 = infected_normal*infected_normal_treatment_rate;  // w[11]
    flux_v8 = infected_moderate*infected_moderate_treatment_rate;  // w[12]
    flux_v9 = infected_weak*infected_weak_treatment_rate;  // w[13]
    flux_v10 = susceptible*susceptible_death_rate;  // w[14]
    flux_v11 = infected_normal*infected_normal_death_rate;  // w[15]
    flux_v12 = infected_moderate*infected_moderate_death_rate;  // w[16]
    flux_v13 = infected_weak*infected_weak_death_rate;  // w[17]
    flux_v14 = treated_normal*treated_normal_death_rate;  // w[18]
    flux_v15 = treated_moderate*treated_moderate_death_rate;  // w[19]
    flux_v16 = treated_weak*treated_weak_death_rate;  // w[20]
}

} // namespace model_Rahman_MBS2016
} // namespace amici
