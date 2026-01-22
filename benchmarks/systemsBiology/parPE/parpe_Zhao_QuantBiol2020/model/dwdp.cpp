#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 26> dwdp_colptrs_Zhao_QuantBiol2020_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 21, 21, 21, 21
};

void dwdp_colptrs_Zhao_QuantBiol2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Zhao_QuantBiol2020_));
}
} // namespace model_Zhao_QuantBiol2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zhao_QuantBiol2020 {

static constexpr std::array<sunindextype, 21> dwdp_rowvals_Zhao_QuantBiol2020_ = {
    1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4, 1, 3, 4
};

void dwdp_rowvals_Zhao_QuantBiol2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Zhao_QuantBiol2020_));
}
} // namespace model_Zhao_QuantBiol2020
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "dwdp.h"

namespace amici {
namespace model_Zhao_QuantBiol2020 {

void dwdp_Zhao_QuantBiol2020(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dR_dR_Stage_I_Wuhan = Trigger_Stage_I*Trigger_Wuhan;  // dwdp[0]
        dgamma_1_dgamma_1_Stage_I_Wuhan = Trigger_Stage_I*Trigger_Wuhan;  // dwdp[1]
        dgamma_2_dgamma_2_Stage_I_Wuhan = Trigger_Stage_I*Trigger_Wuhan;  // dwdp[2]
        dR_dR_Stage_II_Wuhan = Trigger_Stage_II*Trigger_Wuhan;  // dwdp[3]
        dgamma_1_dgamma_1_Stage_II_Wuhan = Trigger_Stage_II*Trigger_Wuhan;  // dwdp[4]
        dgamma_2_dgamma_2_Stage_II_Wuhan = Trigger_Stage_II*Trigger_Wuhan;  // dwdp[5]
        dR_dR_Stage_III_Wuhan = Trigger_Stage_III*Trigger_Wuhan;  // dwdp[6]
        dgamma_1_dgamma_1_Stage_III_Wuhan = Trigger_Stage_III*Trigger_Wuhan;  // dwdp[7]
        dgamma_2_dgamma_2_Stage_III_Wuhan = Trigger_Stage_III*Trigger_Wuhan;  // dwdp[8]
        dR_dR_Stage_I_Hubei = Trigger_Hubei*Trigger_Stage_I;  // dwdp[9]
        dgamma_1_dgamma_1_Stage_I_Hubei = Trigger_Hubei*Trigger_Stage_I;  // dwdp[10]
        dgamma_2_dgamma_2_Stage_I_Hubei = Trigger_Hubei*Trigger_Stage_I;  // dwdp[11]
        dR_dR_Stage_II_Hubei = Trigger_Hubei*Trigger_Stage_II;  // dwdp[12]
        dgamma_1_dgamma_1_Stage_II_Hubei = Trigger_Hubei*Trigger_Stage_II;  // dwdp[13]
        dgamma_2_dgamma_2_Stage_II_Hubei = Trigger_Hubei*Trigger_Stage_II;  // dwdp[14]
        dR_dR_Stage_I_China = Trigger_China*Trigger_Stage_I;  // dwdp[15]
        dgamma_1_dgamma_1_Stage_I_China = Trigger_China*Trigger_Stage_I;  // dwdp[16]
        dgamma_2_dgamma_2_Stage_I_China = Trigger_China*Trigger_Stage_I;  // dwdp[17]
        dR_dR_Stage_II_China = Trigger_China*Trigger_Stage_II;  // dwdp[18]
        dgamma_1_dgamma_1_Stage_II_China = Trigger_China*Trigger_Stage_II;  // dwdp[19]
        dgamma_2_dgamma_2_Stage_II_China = Trigger_China*Trigger_Stage_II;  // dwdp[20]
    }
}

} // namespace model_Zhao_QuantBiol2020
} // namespace amici
