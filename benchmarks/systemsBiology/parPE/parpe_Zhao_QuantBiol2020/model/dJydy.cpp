#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<std::array<sunindextype, 6>, 5> dJydy_colptrs_Zhao_QuantBiol2020_ = {{
    {0, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Zhao_QuantBiol2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Zhao_QuantBiol2020_[index]));
}
} // namespace model_Zhao_QuantBiol2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<std::array<sunindextype, 1>, 5> dJydy_rowvals_Zhao_QuantBiol2020_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Zhao_QuantBiol2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Zhao_QuantBiol2020_[index]));
}
} // namespace model_Zhao_QuantBiol2020
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "k.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"
#include "dJydy.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void dJydy_Zhao_QuantBiol2020(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_confirmed_infected + 1.0*observable_confirmed_infected)/std::pow(sigma_observable_confirmed_infected, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_susceptible + 1.0*observable_susceptible)/std::pow(sigma_observable_susceptible, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_unquarantined_infected + 1.0*observable_unquarantined_infected)/std::pow(sigma_observable_unquarantined_infected, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_quarantined_infected + 1.0*observable_quarantined_infected)/std::pow(sigma_observable_quarantined_infected, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobservable_cumulative_infected + 1.0*observable_cumulative_infected)/std::pow(sigma_observable_cumulative_infected, 2);
            break;
    }
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
