#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<std::array<sunindextype, 3>, 2> dJydy_colptrs_Oliveira_NatCommun2021_ = {{
    {0, 1, 1}, 
    {0, 0, 1}, 
}};

void dJydy_colptrs_Oliveira_NatCommun2021(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Oliveira_NatCommun2021_[index]));
}
} // namespace model_Oliveira_NatCommun2021
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Oliveira_NatCommun2021 {

static constexpr std::array<std::array<sunindextype, 1>, 2> dJydy_rowvals_Oliveira_NatCommun2021_ = {{
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Oliveira_NatCommun2021(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Oliveira_NatCommun2021_[index]));
}
} // namespace model_Oliveira_NatCommun2021
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
namespace model_Oliveira_NatCommun2021 {

void dJydy_Oliveira_NatCommun2021(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*cumulative_deaths - 1.0*mcumulative_deaths)/std::pow(sigma_cumulative_deaths, 2);
            break;
        case 1:
            dJydy[0] = (1.0*cumulative_cases - 1.0*mcumulative_cases)/std::pow(sigma_cumulative_cases, 2);
            break;
    }
}

} // namespace model_Oliveira_NatCommun2021
} // namespace amici
