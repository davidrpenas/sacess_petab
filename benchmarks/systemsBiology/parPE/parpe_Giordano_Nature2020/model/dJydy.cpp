#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<std::array<sunindextype, 8>, 7> dJydy_colptrs_Giordano_Nature2020_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Giordano_Nature2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Giordano_Nature2020_[index]));
}
} // namespace model_Giordano_Nature2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<std::array<sunindextype, 1>, 7> dJydy_rowvals_Giordano_Nature2020_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Giordano_Nature2020(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Giordano_Nature2020_[index]));
}
} // namespace model_Giordano_Nature2020
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
namespace model_Giordano_Nature2020 {

void dJydy_Giordano_Nature2020(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_TotalCases + 1.0*observable_TotalCases)/std::pow(sigma_observable_TotalCases, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_Deaths + 1.0*observable_Deaths)/std::pow(sigma_observable_Deaths, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_Recovered + 1.0*observable_Recovered)/std::pow(sigma_observable_Recovered, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_CurrentCases + 1.0*observable_CurrentCases)/std::pow(sigma_observable_CurrentCases, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobservable_DiagHome + 1.0*observable_DiagHome)/std::pow(sigma_observable_DiagHome, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobservable_Hospit + 1.0*observable_Hospit)/std::pow(sigma_observable_Hospit, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mobservable_ICU + 1.0*observable_ICU)/std::pow(sigma_observable_ICU, 2);
            break;
    }
}

} // namespace model_Giordano_Nature2020
} // namespace amici
