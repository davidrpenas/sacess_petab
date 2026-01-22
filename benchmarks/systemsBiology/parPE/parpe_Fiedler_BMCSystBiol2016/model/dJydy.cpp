#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<std::array<sunindextype, 13>, 12> dJydy_colptrs_Fiedler_BMCSystBiol2016_ = {{
    {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}, 
}};

void dJydy_colptrs_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexptrs(gsl::make_span(dJydy_colptrs_Fiedler_BMCSystBiol2016_[index]));
}
} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<std::array<sunindextype, 1>, 12> dJydy_rowvals_Fiedler_BMCSystBiol2016_ = {{
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

void dJydy_rowvals_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dJydy, int index){
    dJydy.set_indexvals(gsl::make_span(dJydy_rowvals_Fiedler_BMCSystBiol2016_[index]));
}
} // namespace model_Fiedler_BMCSystBiol2016
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
namespace model_Fiedler_BMCSystBiol2016 {

void dJydy_Fiedler_BMCSystBiol2016(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mpErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1 + 1.0*pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1)/std::pow(sigma_pErk__s_pErk_20140430_gel1__s_pErk_20140430_gel1_sigma_pErk__model1_data1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mpErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2 + 1.0*pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2)/std::pow(sigma_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data2, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mpErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3 + 1.0*pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3)/std::pow(sigma_pErk__s_pErk_20140430_gel2__s_pErk_20140430_gel2_sigma_pErk__model1_data3, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mpErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1 + 1.0*pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1)/std::pow(sigma_pErk__s_pErk_20140505_gel1__s_pErk_20140505_gel1_sigma_pErk__model1_data1, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mpErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2 + 1.0*pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2)/std::pow(sigma_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data2, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mpErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3 + 1.0*pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3)/std::pow(sigma_pErk__s_pErk_20140505_gel2__s_pErk_20140505_gel2_sigma_pErk__model1_data3, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mpMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1 + 1.0*pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1)/std::pow(sigma_pMek__s_pMek_20140430_gel1__s_pMek_20140430_gel1_sigma_pMek__model1_data1, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mpMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2 + 1.0*pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2)/std::pow(sigma_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data2, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mpMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3 + 1.0*pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3)/std::pow(sigma_pMek__s_pMek_20140430_gel2__s_pMek_20140430_gel2_sigma_pMek__model1_data3, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*mpMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1 + 1.0*pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1)/std::pow(sigma_pMek__s_pMek_20140505_gel1__s_pMek_20140505_gel1_sigma_pMek__model1_data1, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*mpMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2 + 1.0*pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2)/std::pow(sigma_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data2, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*mpMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3 + 1.0*pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3)/std::pow(sigma_pMek__s_pMek_20140505_gel2__s_pMek_20140505_gel2_sigma_pMek__model1_data3, 2);
            break;
    }
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
