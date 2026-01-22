#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<sunindextype, 19> dwdp_colptrs_Armistead_CellDeathDis2024_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10
};

void dwdp_colptrs_Armistead_CellDeathDis2024(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Armistead_CellDeathDis2024_));
}
} // namespace model_Armistead_CellDeathDis2024
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<sunindextype, 10> dwdp_rowvals_Armistead_CellDeathDis2024_ = {
    0, 1, 2, 3, 5, 6, 7, 4, 0, 5
};

void dwdp_rowvals_Armistead_CellDeathDis2024(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Armistead_CellDeathDis2024_));
}
} // namespace model_Armistead_CellDeathDis2024
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
namespace model_Armistead_CellDeathDis2024 {

void dwdp_Armistead_CellDeathDis2024(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_R1_dk00 = S_on*alpha_cer + 1;  // dwdp[0]
        dflux_R1_dalpha_cer = S_on*k00;  // dwdp[8]
    }

    // dynamic expressions
    dflux_R2_dk0 = Sphinga;  // dwdp[1]
    dflux_R3_dk1 = Cer;  // dwdp[2]
    dflux_R4_dk2 = Sphingo;  // dwdp[3]
    dflux_R6_dk3 = Sphingo*(-S_on*alpha_hai1a + 1);  // dwdp[4]
    dflux_R7_dk4 = S1P;  // dwdp[5]
    dflux_R8_dk5 = S1P;  // dwdp[6]
    dflux_R5_dk_d = Cer;  // dwdp[7]
    dflux_R6_dalpha_hai1a = -S_on*Sphingo*k3;  // dwdp[9]
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
