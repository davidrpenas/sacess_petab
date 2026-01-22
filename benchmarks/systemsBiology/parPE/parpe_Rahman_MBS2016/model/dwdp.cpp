#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 11> dwdp_colptrs_Rahman_MBS2016_ = {
    0, 1, 5, 6, 7, 8, 9, 10, 11, 12, 12
};

void dwdp_colptrs_Rahman_MBS2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Rahman_MBS2016_));
}
} // namespace model_Rahman_MBS2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 12> dwdp_rowvals_Rahman_MBS2016_ = {
    0, 0, 1, 3, 4, 1, 13, 7, 8, 10, 9, 4
};

void dwdp_rowvals_Rahman_MBS2016(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Rahman_MBS2016_));
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
#include "dwdp.h"

namespace amici {
namespace model_Rahman_MBS2016 {

void dwdp_Rahman_MBS2016(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dinfected_normal_transmission_rate_dinfected_normal_transmission_rate_relative = infected_moderate_transmission_rate;  // dwdp[0]
        dinfected_normal_transmission_rate_dinfected_moderate_transmission_rate = infected_normal_transmission_rate_relative;  // dwdp[1]
        dinfected_weak_transmission_rate_dinfected_moderate_transmission_rate = infected_weak_transmission_rate_relative;  // dwdp[2]
        dtreated_transmission_rate_dinfected_moderate_transmission_rate = 0.040000000000000001;  // dwdp[3]
        dinfected_weak_transmission_rate_dinfected_weak_transmission_rate_relative = infected_moderate_transmission_rate;  // dwdp[5]
    }

    // dynamic expressions
    dforce_of_infection_dinfected_moderate_transmission_rate = infected_moderate*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);  // dwdp[4]
    dflux_v9_dinfected_weak_treatment_rate = infected_weak;  // dwdp[6]
    dflux_v3_dinfected_normal_worsen_rate = infected_normal;  // dwdp[7]
    dflux_v4_dinfected_moderate_worsen_rate = infected_moderate;  // dwdp[8]
    dflux_v6_dtreated_moderate_improve_rate = treated_moderate;  // dwdp[9]
    dflux_v5_dtreated_weak_improve_rate = treated_weak;  // dwdp[10]
    dforce_of_infection_dbehavioural_change_rate = (infected_moderate*infected_moderate_transmission_rate + infected_normal*infected_normal_transmission_rate + infected_weak*infected_weak_transmission_rate + treated_transmission_rate*(treated_moderate + treated_normal + treated_weak))*(-infected_moderate - infected_normal - infected_weak - treated_moderate - treated_normal - treated_weak)*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);  // dwdp[11]
}

} // namespace model_Rahman_MBS2016
} // namespace amici
