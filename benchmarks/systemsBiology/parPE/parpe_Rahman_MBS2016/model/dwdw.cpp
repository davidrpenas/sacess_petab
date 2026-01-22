#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 22> dwdw_colptrs_Rahman_MBS2016_ = {
    0, 1, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

void dwdw_colptrs_Rahman_MBS2016(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Rahman_MBS2016_));
}
} // namespace model_Rahman_MBS2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<sunindextype, 4> dwdw_rowvals_Rahman_MBS2016_ = {
    4, 4, 4, 6
};

void dwdw_rowvals_Rahman_MBS2016(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Rahman_MBS2016_));
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
#include "dwdw.h"

namespace amici {
namespace model_Rahman_MBS2016 {

void dwdw_Rahman_MBS2016(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dforce_of_infection_dinfected_normal_transmission_rate = infected_normal*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);  // dwdw[0]
    dforce_of_infection_dinfected_weak_transmission_rate = infected_weak*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);  // dwdw[1]
    dforce_of_infection_dtreated_transmission_rate = (treated_moderate + treated_normal + treated_weak)*std::exp(-behavioural_change_rate*(infected_moderate + infected_normal + infected_weak + treated_moderate + treated_normal + treated_weak))/(infected_moderate + infected_normal + infected_weak + susceptible + treated_moderate + treated_normal + treated_weak);  // dwdw[2]
    dflux_v2_dforce_of_infection = susceptible;  // dwdw[3]
}

} // namespace model_Rahman_MBS2016
} // namespace amici
