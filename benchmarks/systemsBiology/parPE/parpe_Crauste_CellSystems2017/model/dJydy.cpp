#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<std::array<sunindextype, 22>, 21> dJydy_colptrs_Crauste_CellSystems2017_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Crauste_CellSystems2017(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Crauste_CellSystems2017_[index]));
}
} // namespace model_Crauste_CellSystems2017
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Crauste_CellSystems2017 {

static constexpr std::array<std::array<sunindextype, 1>, 21> dJydy_rowvals_Crauste_CellSystems2017_ = {{
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
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Crauste_CellSystems2017(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Crauste_CellSystems2017_[index]));
}
} // namespace model_Crauste_CellSystems2017
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
namespace model_Crauste_CellSystems2017 {

void dJydy_Crauste_CellSystems2017(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_EarlyEffector__266__model1_data1 + 1.0*observable_EarlyEffector__266__model1_data1)/std::pow(sigma_observable_EarlyEffector__266__model1_data1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_EarlyEffector__526__model1_data1 + 1.0*observable_EarlyEffector__526__model1_data1)/std::pow(sigma_observable_EarlyEffector__526__model1_data1, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_EarlyEffector__1540__model1_data1 + 1.0*observable_EarlyEffector__1540__model1_data1)/std::pow(sigma_observable_EarlyEffector__1540__model1_data1, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_EarlyEffector__2924__model1_data1 + 1.0*observable_EarlyEffector__2924__model1_data1)/std::pow(sigma_observable_EarlyEffector__2924__model1_data1, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobservable_EarlyEffector__3048__model1_data1 + 1.0*observable_EarlyEffector__3048__model1_data1)/std::pow(sigma_observable_EarlyEffector__3048__model1_data1, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobservable_EarlyEffector__19722__model1_data1 + 1.0*observable_EarlyEffector__19722__model1_data1)/std::pow(sigma_observable_EarlyEffector__19722__model1_data1, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mobservable_EarlyEffector__23855__model1_data1 + 1.0*observable_EarlyEffector__23855__model1_data1)/std::pow(sigma_observable_EarlyEffector__23855__model1_data1, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mobservable_LateEffector__182__model1_data1 + 1.0*observable_LateEffector__182__model1_data1)/std::pow(sigma_observable_LateEffector__182__model1_data1, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mobservable_LateEffector__1976__model1_data1 + 1.0*observable_LateEffector__1976__model1_data1)/std::pow(sigma_observable_LateEffector__1976__model1_data1, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*mobservable_LateEffector__2705__model1_data1 + 1.0*observable_LateEffector__2705__model1_data1)/std::pow(sigma_observable_LateEffector__2705__model1_data1, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*mobservable_LateEffector__4105__model1_data1 + 1.0*observable_LateEffector__4105__model1_data1)/std::pow(sigma_observable_LateEffector__4105__model1_data1, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*mobservable_LateEffector__6841__model1_data1 + 1.0*observable_LateEffector__6841__model1_data1)/std::pow(sigma_observable_LateEffector__6841__model1_data1, 2);
            break;
        case 12:
            dJydy[0] = (-1.0*mobservable_LateEffector__7266__model1_data1 + 1.0*observable_LateEffector__7266__model1_data1)/std::pow(sigma_observable_LateEffector__7266__model1_data1, 2);
            break;
        case 13:
            dJydy[0] = (-1.0*mobservable_LateEffector__23601__model1_data1 + 1.0*observable_LateEffector__23601__model1_data1)/std::pow(sigma_observable_LateEffector__23601__model1_data1, 2);
            break;
        case 14:
            dJydy[0] = (-1.0*mobservable_Memory__80__model1_data1 + 1.0*observable_Memory__80__model1_data1)/std::pow(sigma_observable_Memory__80__model1_data1, 2);
            break;
        case 15:
            dJydy[0] = (-1.0*mobservable_Memory__4206__model1_data1 + 1.0*observable_Memory__4206__model1_data1)/std::pow(sigma_observable_Memory__4206__model1_data1, 2);
            break;
        case 16:
            dJydy[0] = (-1.0*mobservable_Memory__5745__model1_data1 + 1.0*observable_Memory__5745__model1_data1)/std::pow(sigma_observable_Memory__5745__model1_data1, 2);
            break;
        case 17:
            dJydy[0] = (-1.0*mobservable_Memory__7909__model1_data1 + 1.0*observable_Memory__7909__model1_data1)/std::pow(sigma_observable_Memory__7909__model1_data1, 2);
            break;
        case 18:
            dJydy[0] = (-1.0*mobservable_Memory__15987__model1_data1 + 1.0*observable_Memory__15987__model1_data1)/std::pow(sigma_observable_Memory__15987__model1_data1, 2);
            break;
        case 19:
            dJydy[0] = (-1.0*mobservable_Naive__64__model1_data1 + 1.0*observable_Naive__64__model1_data1)/std::pow(sigma_observable_Naive__64__model1_data1, 2);
            break;
        case 20:
            dJydy[0] = (-1.0*mobservable_Naive__290__model1_data1 + 1.0*observable_Naive__290__model1_data1)/std::pow(sigma_observable_Naive__290__model1_data1, 2);
            break;
    }
}

} // namespace model_Crauste_CellSystems2017
} // namespace amici
