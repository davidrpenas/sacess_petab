#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<std::array<sunindextype, 4>, 3> dJydy_colptrs_Boehm_JProteomeRes2014_ = {{
    {0, 1, 1, 1}, 
    {0, 0, 1, 1}, 
    {0, 0, 0, 1}, 
}};

void dJydy_colptrs_Boehm_JProteomeRes2014(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Boehm_JProteomeRes2014_[index]));
}
} // namespace model_Boehm_JProteomeRes2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<std::array<sunindextype, 1>, 3> dJydy_rowvals_Boehm_JProteomeRes2014_ = {{
    {0}, 
    {0}, 
    {0}, 
}};

void dJydy_rowvals_Boehm_JProteomeRes2014(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Boehm_JProteomeRes2014_[index]));
}
} // namespace model_Boehm_JProteomeRes2014
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
namespace model_Boehm_JProteomeRes2014 {

void dJydy_Boehm_JProteomeRes2014(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mpSTAT5A_rel + 1.0*pSTAT5A_rel)/std::pow(sigma_pSTAT5A_rel, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mpSTAT5B_rel + 1.0*pSTAT5B_rel)/std::pow(sigma_pSTAT5B_rel, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mrSTAT5A_rel + 1.0*rSTAT5A_rel)/std::pow(sigma_rSTAT5A_rel, 2);
            break;
    }
}

} // namespace model_Boehm_JProteomeRes2014
} // namespace amici
