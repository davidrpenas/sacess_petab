#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 8> dwdx_colptrs_Rahman_MBS2016_ = {
    0, 4, 9, 14, 18, 21, 25, 29
};

void dwdx_colptrs_Rahman_MBS2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Rahman_MBS2016_));
}
} // namespace model_Rahman_MBS2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 29> dwdx_rowvals_Rahman_MBS2016_ = {
    2, 4, 6, 14, 2, 4, 7, 11, 15, 2, 4, 8, 12, 16, 2, 4, 13, 17, 2, 4, 18, 2, 4, 10, 19, 2, 4, 9, 20
};

void dwdx_rowvals_Rahman_MBS2016(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Rahman_MBS2016_));
}
} // namespace model_Rahman_MBS2016
} // namespace amici




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

void dwdx_Rahman_MBS2016(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v10_dsusceptible = susceptible_death_rate;  // dwdx[3]
        dflux_v3_dinfected_normal = infected_normal_worsen_rate;  // dwdx[6]
        dflux_v7_dinfected_normal = infected_normal_treatment_rate;  // dwdx[7]
        dflux_v11_dinfected_normal = infected_normal_death_rate;  // dwdx[8]
        dflux_v4_dinfected_moderate = infected_moderate_worsen_rate;  // dwdx[11]
        dflux_v8_dinfected_moderate = infected_moderate_treatment_rate;  // dwdx[12]
        dflux_v12_dinfected_moderate = infected_moderate_death_rate;  // dwdx[13]
        dflux_v9_dinfected_weak = infected_weak_treatment_rate;  // dwdx[16]
        dflux_v13_dinfected_weak = infected_weak_death_rate;  // dwdx[17]
        dflux_v14_dtreated_normal = treated_normal_death_rate;  // dwdx[20]
        dflux_v6_dtreated_moderate = treated_moderate_improve_rate;  // dwdx[23]
        dflux_v15_dtreated_moderate = treated_moderate_death_rate;  // dwdx[24]
        dflux_v5_dtreated_weak = treated_weak_improve_rate;  // dwdx[27]
        dflux_v16_dtreated_weak = treated_weak_death_rate;  // dwdx[28]
    }

    // dynamic expressions
    dprevalence_dsusceptible = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2) - 100/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);  // dwdx[0]
    dforce_of_infection_dsusceptible = (-infected_moderate*infected_moderate_transmission_rate - infected_normal*infected_normal_transmission_rate - infected_weak*infected_weak_transmission_rate - treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[1]
    dflux_v2_dsusceptible = force_of_infection;  // dwdx[2]
    dprevalence_dinfected_normal = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[4]
    dforce_of_infection_dinfected_normal = -behavioural_change_rate*(infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + infected_normal_transmission_rate*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + (-infected_moderate*infected_moderate_transmission_rate - infected_normal*infected_normal_transmission_rate - infected_weak*infected_weak_transmission_rate - treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[5]
    dprevalence_dinfected_moderate = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[9]
    dforce_of_infection_dinfected_moderate = -behavioural_change_rate*(infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + infected_moderate_transmission_rate*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + (-infected_moderate*infected_moderate_transmission_rate - infected_normal*infected_normal_transmission_rate - infected_weak*infected_weak_transmission_rate - treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[10]
    dprevalence_dinfected_weak = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[14]
    dforce_of_infection_dinfected_weak = -behavioural_change_rate*(infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + infected_weak_transmission_rate*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + (-infected_moderate*infected_moderate_transmission_rate - infected_normal*infected_normal_transmission_rate - infected_weak*infected_weak_transmission_rate - treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[15]
    dprevalence_dtreated_normal = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[18]
    dforce_of_infection_dtreated_normal = -behavioural_change_rate*(infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + treated_transmission_rate*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + (-infected_moderate*infected_moderate_transmission_rate - infected_normal*infected_normal_transmission_rate - infected_weak*infected_weak_transmission_rate - treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[19]
    dprevalence_dtreated_moderate = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[21]
    dforce_of_infection_dtreated_moderate = -behavioural_change_rate*(infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + treated_transmission_rate*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + (-infected_moderate*infected_moderate_transmission_rate - infected_normal*infected_normal_transmission_rate - infected_weak*infected_weak_transmission_rate - treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[22]
    dprevalence_dtreated_weak = 100*susceptible/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[25]
    dforce_of_infection_dtreated_weak = -behavioural_change_rate*(infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + treated_transmission_rate*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak) + (-infected_moderate*infected_moderate_transmission_rate - infected_normal*infected_normal_transmission_rate - infected_weak*infected_weak_transmission_rate - treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/std::pow(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak, 2);  // dwdx[26]
}

} // namespace model_Rahman_MBS2016
} // namespace amici
