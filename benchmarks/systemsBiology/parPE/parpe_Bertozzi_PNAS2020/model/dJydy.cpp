#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<std::array<sunindextype, 3>, 2> dJydy_colptrs_Bertozzi_PNAS2020_ = {{
    {0, 1, 1}, 
    {0, 0, 1}, 
}};

void dJydy_colptrs_Bertozzi_PNAS2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Bertozzi_PNAS2020_[index]));
}
} // namespace model_Bertozzi_PNAS2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Bertozzi_PNAS2020 {

static constexpr std::array<std::array<sunindextype, 1>, 2> dJydy_rowvals_Bertozzi_PNAS2020_ = {{
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Bertozzi_PNAS2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Bertozzi_PNAS2020_[index]));
}
} // namespace model_Bertozzi_PNAS2020
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
namespace model_Bertozzi_PNAS2020 {

void dJydy_Bertozzi_PNAS2020(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*my_I_NY + 1.0*y_I_NY)/std::pow(sigma_y_I_NY, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*my_I_CA + 1.0*y_I_CA)/std::pow(sigma_y_I_CA, 2);
            break;
    }
}

} // namespace model_Bertozzi_PNAS2020
} // namespace amici
