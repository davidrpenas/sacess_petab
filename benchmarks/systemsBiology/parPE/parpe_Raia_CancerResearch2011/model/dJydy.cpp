#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<std::array<sunindextype, 9>, 8> dJydy_colptrs_Raia_CancerResearch2011_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Raia_CancerResearch2011(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Raia_CancerResearch2011_[index]));
}
} // namespace model_Raia_CancerResearch2011
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Raia_CancerResearch2011 {

static constexpr std::array<std::array<sunindextype, 1>, 8> dJydy_rowvals_Raia_CancerResearch2011_ = {{
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Raia_CancerResearch2011(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Raia_CancerResearch2011_[index]));
}
} // namespace model_Raia_CancerResearch2011
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
namespace model_Raia_CancerResearch2011 {

void dJydy_Raia_CancerResearch2011(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_CD274mRNA + 1.0*observable_CD274mRNA)/std::pow(sigma_observable_CD274mRNA, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_IL13_cell + 1.0*observable_IL13_cell)/std::pow(sigma_observable_IL13_cell, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_RecSurf + 1.0*observable_RecSurf)/std::pow(sigma_observable_RecSurf, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_SOCS3 + 1.0*observable_SOCS3)/std::pow(sigma_observable_SOCS3, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobservable_SOCS3mRNA + 1.0*observable_SOCS3mRNA)/std::pow(sigma_observable_SOCS3mRNA, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobservable_pIL4Ra + 1.0*observable_pIL4Ra)/std::pow(sigma_observable_pIL4Ra, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mobservable_pJAK2 + 1.0*observable_pJAK2)/std::pow(sigma_observable_pJAK2, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mobservable_pSTAT5 + 1.0*observable_pSTAT5)/std::pow(sigma_observable_pSTAT5, 2);
            break;
    }
}

} // namespace model_Raia_CancerResearch2011
} // namespace amici
