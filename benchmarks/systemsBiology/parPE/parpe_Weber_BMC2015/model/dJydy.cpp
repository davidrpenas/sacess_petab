#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<std::array<sunindextype, 9>, 8> dJydy_colptrs_Weber_BMC2015_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Weber_BMC2015(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Weber_BMC2015_[index]));
}
} // namespace model_Weber_BMC2015
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<std::array<sunindextype, 1>, 8> dJydy_rowvals_Weber_BMC2015_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Weber_BMC2015(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Weber_BMC2015_[index]));
}
} // namespace model_Weber_BMC2015
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
namespace model_Weber_BMC2015 {

void dJydy_Weber_BMC2015(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*myCERTpRN24 + 1.0*yCERTpRN24)/std::pow(sigma_yCERTpRN24, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*myCERTt + 1.0*yCERTt)/std::pow(sigma_yCERTt, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*myPI4K3BpRN24 + 1.0*yPI4K3BpRN24)/std::pow(sigma_yPI4K3BpRN24, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*myPI4K3Bt + 1.0*yPI4K3Bt)/std::pow(sigma_yPI4K3Bt, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*myPKDpN0 + 1.0*yPKDpN0)/std::pow(sigma_yPKDpN0, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*myPKDpN24 + 1.0*yPKDpN24)/std::pow(sigma_yPKDpN24, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*myPKDpN25 + 1.0*yPKDpN25)/std::pow(sigma_yPKDpN25, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*myPKDt + 1.0*yPKDt)/std::pow(sigma_yPKDt, 2);
            break;
    }
}

} // namespace model_Weber_BMC2015
} // namespace amici
