#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<std::array<sunindextype, 4>, 3> dJydy_colptrs_Brannmark_JBC2010_ = {{
    {0, 1, 1, 1}, 
    {0, 0, 1, 1}, 
    {0, 0, 0, 1}, 
}};

void dJydy_colptrs_Brannmark_JBC2010(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Brannmark_JBC2010_[index]));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<std::array<sunindextype, 1>, 3> dJydy_rowvals_Brannmark_JBC2010_ = {{
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Brannmark_JBC2010_[index]));
}
} // namespace model_Brannmark_JBC2010
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
namespace model_Brannmark_JBC2010 {

void dJydy_Brannmark_JBC2010(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (1.0*IR1_P - 1.0*mIR1_P)/std::pow(sigma_IR1_P, 2);
            break;
        case 1:
            dJydy[0] = (1.0*IRS1_P - 1.0*mIRS1_P)/std::pow(sigma_IRS1_P, 2);
            break;
        case 2:
            dJydy[0] = (1.0*IRS1_P_DosR - 1.0*mIRS1_P_DosR)/std::pow(sigma_IRS1_P_DosR, 2);
            break;
    }
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
