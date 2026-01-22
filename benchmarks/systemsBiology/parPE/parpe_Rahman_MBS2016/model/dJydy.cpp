#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<std::array<sunindextype, 2>, 1> dJydy_colptrs_Rahman_MBS2016_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_Rahman_MBS2016(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Rahman_MBS2016_[index]));
}
} // namespace model_Rahman_MBS2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Rahman_MBS2016 {

static constexpr std::array<std::array<sunindextype, 1>, 1> dJydy_rowvals_Rahman_MBS2016_ = {{
    {0}, 
}};

void dJydy_rowvals_Rahman_MBS2016(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Rahman_MBS2016_[index]));
}
} // namespace model_Rahman_MBS2016
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
namespace model_Rahman_MBS2016 {

void dJydy_Rahman_MBS2016(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_prevalence + 1.0*observable_prevalence)/std::pow(sigma_observable_prevalence, 2);
            break;
    }
}

} // namespace model_Rahman_MBS2016
} // namespace amici
