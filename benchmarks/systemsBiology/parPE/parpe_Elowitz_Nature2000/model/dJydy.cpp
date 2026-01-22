#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<std::array<sunindextype, 2>, 1> dJydy_colptrs_Elowitz_Nature2000_ = {{
    {0, 1}, 
}};

void dJydy_colptrs_Elowitz_Nature2000(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Elowitz_Nature2000_[index]));
}
} // namespace model_Elowitz_Nature2000
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Elowitz_Nature2000 {

static constexpr std::array<std::array<sunindextype, 1>, 1> dJydy_rowvals_Elowitz_Nature2000_ = {{
    {0}, 
}};

void dJydy_rowvals_Elowitz_Nature2000(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Elowitz_Nature2000_[index]));
}
} // namespace model_Elowitz_Nature2000
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"
#include "dJydy.h"

namespace amici {
namespace model_Elowitz_Nature2000 {

void dJydy_Elowitz_Nature2000(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*std::log(fluorescence)/M_LN10 - 1.0*std::log(mfluorescence)/M_LN10)/(fluorescence*std::pow(sigma_fluorescence, 2)*M_LN10);
            break;
    }
}

} // namespace model_Elowitz_Nature2000
} // namespace amici
