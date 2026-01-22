#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<std::array<sunindextype, 17>, 16> dJydy_colptrs_Blasi_CellSystems2016_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Blasi_CellSystems2016(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Blasi_CellSystems2016_[index]));
}
} // namespace model_Blasi_CellSystems2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Blasi_CellSystems2016 {

static constexpr std::array<std::array<sunindextype, 1>, 16> dJydy_rowvals_Blasi_CellSystems2016_ = {{
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

void dJydy_rowvals_Blasi_CellSystems2016(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Blasi_CellSystems2016_[index]));
}
} // namespace model_Blasi_CellSystems2016
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
namespace model_Blasi_CellSystems2016 {

void dJydy_Blasi_CellSystems2016(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*std::log(mobservable_0ac) + 1.0*std::log(observable_0ac))/(observable_0ac*std::pow(sigma_observable_0ac, 2));
            break;
        case 1:
            dJydy[0] = (-1.0*std::log(mobservable_4ac) + 1.0*std::log(observable_4ac))/(observable_4ac*std::pow(sigma_observable_4ac, 2));
            break;
        case 2:
            dJydy[0] = (-1.0*std::log(mobservable_k12) + 1.0*std::log(observable_k12))/(observable_k12*std::pow(sigma_observable_k12, 2));
            break;
        case 3:
            dJydy[0] = (-1.0*std::log(mobservable_k12k16) + 1.0*std::log(observable_k12k16))/(observable_k12k16*std::pow(sigma_observable_k12k16, 2));
            break;
        case 4:
            dJydy[0] = (-1.0*std::log(mobservable_k16) + 1.0*std::log(observable_k16))/(observable_k16*std::pow(sigma_observable_k16, 2));
            break;
        case 5:
            dJydy[0] = (-1.0*std::log(mobservable_k5) + 1.0*std::log(observable_k5))/(observable_k5*std::pow(sigma_observable_k5, 2));
            break;
        case 6:
            dJydy[0] = (-1.0*std::log(mobservable_k5k12) + 1.0*std::log(observable_k5k12))/(observable_k5k12*std::pow(sigma_observable_k5k12, 2));
            break;
        case 7:
            dJydy[0] = (-1.0*std::log(mobservable_k5k12k16) + 1.0*std::log(observable_k5k12k16))/(observable_k5k12k16*std::pow(sigma_observable_k5k12k16, 2));
            break;
        case 8:
            dJydy[0] = (-1.0*std::log(mobservable_k5k16) + 1.0*std::log(observable_k5k16))/(observable_k5k16*std::pow(sigma_observable_k5k16, 2));
            break;
        case 9:
            dJydy[0] = (-1.0*std::log(mobservable_k5k8) + 1.0*std::log(observable_k5k8))/(observable_k5k8*std::pow(sigma_observable_k5k8, 2));
            break;
        case 10:
            dJydy[0] = (-1.0*std::log(mobservable_k5k8k12) + 1.0*std::log(observable_k5k8k12))/(observable_k5k8k12*std::pow(sigma_observable_k5k8k12, 2));
            break;
        case 11:
            dJydy[0] = (-1.0*std::log(mobservable_k5k8k16) + 1.0*std::log(observable_k5k8k16))/(observable_k5k8k16*std::pow(sigma_observable_k5k8k16, 2));
            break;
        case 12:
            dJydy[0] = (-1.0*std::log(mobservable_k8) + 1.0*std::log(observable_k8))/(observable_k8*std::pow(sigma_observable_k8, 2));
            break;
        case 13:
            dJydy[0] = (-1.0*std::log(mobservable_k8k12) + 1.0*std::log(observable_k8k12))/(observable_k8k12*std::pow(sigma_observable_k8k12, 2));
            break;
        case 14:
            dJydy[0] = (-1.0*std::log(mobservable_k8k12k16) + 1.0*std::log(observable_k8k12k16))/(observable_k8k12k16*std::pow(sigma_observable_k8k12k16, 2));
            break;
        case 15:
            dJydy[0] = (-1.0*std::log(mobservable_k8k16) + 1.0*std::log(observable_k8k16))/(observable_k8k16*std::pow(sigma_observable_k8k16, 2));
            break;
    }
}

} // namespace model_Blasi_CellSystems2016
} // namespace amici
