#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<std::array<sunindextype, 16>, 15> dJydy_colptrs_Zheng_PNAS2012_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Zheng_PNAS2012(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Zheng_PNAS2012_[index]));
}
} // namespace model_Zheng_PNAS2012
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<std::array<sunindextype, 1>, 15> dJydy_rowvals_Zheng_PNAS2012_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Zheng_PNAS2012(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Zheng_PNAS2012_[index]));
}
} // namespace model_Zheng_PNAS2012
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
namespace model_Zheng_PNAS2012 {

void dJydy_Zheng_PNAS2012(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_K27me0K36me0 + 1.0*observable_K27me0K36me0)/std::pow(sigma_observable_K27me0K36me0, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_K27me0K36me1 + 1.0*observable_K27me0K36me1)/std::pow(sigma_observable_K27me0K36me1, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_K27me0K36me2 + 1.0*observable_K27me0K36me2)/std::pow(sigma_observable_K27me0K36me2, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_K27me0K36me3 + 1.0*observable_K27me0K36me3)/std::pow(sigma_observable_K27me0K36me3, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobservable_K27me1K36me0 + 1.0*observable_K27me1K36me0)/std::pow(sigma_observable_K27me1K36me0, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobservable_K27me1K36me1 + 1.0*observable_K27me1K36me1)/std::pow(sigma_observable_K27me1K36me1, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mobservable_K27me1K36me2 + 1.0*observable_K27me1K36me2)/std::pow(sigma_observable_K27me1K36me2, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mobservable_K27me1K36me3 + 1.0*observable_K27me1K36me3)/std::pow(sigma_observable_K27me1K36me3, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mobservable_K27me2K36me0 + 1.0*observable_K27me2K36me0)/std::pow(sigma_observable_K27me2K36me0, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*mobservable_K27me2K36me1 + 1.0*observable_K27me2K36me1)/std::pow(sigma_observable_K27me2K36me1, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*mobservable_K27me2K36me2 + 1.0*observable_K27me2K36me2)/std::pow(sigma_observable_K27me2K36me2, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*mobservable_K27me2K36me3 + 1.0*observable_K27me2K36me3)/std::pow(sigma_observable_K27me2K36me3, 2);
            break;
        case 12:
            dJydy[0] = (-1.0*mobservable_K27me3K36me0 + 1.0*observable_K27me3K36me0)/std::pow(sigma_observable_K27me3K36me0, 2);
            break;
        case 13:
            dJydy[0] = (-1.0*mobservable_K27me3K36me1 + 1.0*observable_K27me3K36me1)/std::pow(sigma_observable_K27me3K36me1, 2);
            break;
        case 14:
            dJydy[0] = (-1.0*mobservable_K27me3K36me2 + 1.0*observable_K27me3K36me2)/std::pow(sigma_observable_K27me3K36me2, 2);
            break;
    }
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
