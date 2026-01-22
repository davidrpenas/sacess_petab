#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 40> dwdp_colptrs_Weber_BMC2015_ = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28
};

void dwdp_colptrs_Weber_BMC2015(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Weber_BMC2015_));
}
} // namespace model_Weber_BMC2015
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 28> dwdp_rowvals_Weber_BMC2015_ = {
    10, 11, 16, 17, 23, 24, 25, 5, 13, 5, 18, 20, 5, 6, 7, 12, 13, 5, 18, 19, 20, 15, 22, 6, 7, 8, 14, 21
};

void dwdp_rowvals_Weber_BMC2015(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Weber_BMC2015_));
}
} // namespace model_Weber_BMC2015
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"
#include "h.h"
#include "w.h"
#include "dwdp.h"

namespace amici {
namespace model_Weber_BMC2015 {

void dwdp_Weber_BMC2015(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v4_v_3_ds12 = 1.0;  // dwdp[25]
        dflux_v10_v_9_ds21 = 1.0;  // dwdp[26]
        dflux_v17_v_16_ds31 = 1.0;  // dwdp[27]
    }

    // dynamic expressions
    dflux_v6_v_5_da11 = 1.0*PKD;  // dwdp[0]
    dflux_v7_v_6_da12 = 1.0*PKDDAGa;  // dwdp[1]
    dflux_v12_v_11_da21 = 1.0*PI4K3B;  // dwdp[2]
    dflux_v13_v_12_da22 = 1.0*PI4K3Ba;  // dwdp[3]
    dflux_v19_v_18_da31 = 1.0*CERTERa;  // dwdp[4]
    dflux_v20_v_19_da32 = 1.0*CERT;  // dwdp[5]
    dflux_v21_v_20_da33 = 1.0*CERTTGNa;  // dwdp[6]
    dflux_v1_v_0_dm11 = -1.0*CERTERa*PI4K3Ba*PKD*p11*p31/((PI4K3Ba + m31)*std::pow(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11, 2));  // dwdp[7]
    dflux_v9_v_8_dm22 = -1.0*PI4K3B*PKDDAGa*p22/std::pow(PKDDAGa + m22, 2);  // dwdp[8]
    dflux_v1_v_0_dm31 = 1.0*std::pow(CERTERa, 2)*std::pow(PI4K3Ba, 2)*PKD*p11*std::pow(p31, 2)/(std::pow(PI4K3Ba + m31, 3)*std::pow(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11, 2)) - 1.0*CERTERa*PI4K3Ba*PKD*p11*p31/(std::pow(PI4K3Ba + m31, 2)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11));  // dwdp[9]
    dflux_v14_v_13_dm31 = -1.0*CERTERa*PI4K3Ba*p31/std::pow(PI4K3Ba + m31, 2);  // dwdp[10]
    dflux_v16_v_15_dm33 = -1.0*CERTTGNa*PKDDAGa*p33/std::pow(PKDDAGa + m33, 2);  // dwdp[11]
    dflux_v1_v_0_dp11 = 1.0*CERTERa*PI4K3Ba*PKD*p31/((PI4K3Ba + m31)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11));  // dwdp[12]
    dflux_v2_v_1_dp12 = 1.0*PKD*(pu5*u5 + 1);  // dwdp[13]
    dflux_v3_v_2_dp13 = 1.0*PKDDAGa*(pu6*u6 + 1);  // dwdp[14]
    dflux_v8_v_7_dp21 = 1.0*PI4K3Ba;  // dwdp[15]
    dflux_v9_v_8_dp22 = 1.0*PI4K3B*PKDDAGa/(PKDDAGa + m22);  // dwdp[16]
    dflux_v1_v_0_dp31 = -1.0*std::pow(CERTERa, 2)*std::pow(PI4K3Ba, 2)*PKD*p11*p31/(std::pow(PI4K3Ba + m31, 2)*std::pow(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11, 2)) + 1.0*CERTERa*PI4K3Ba*PKD*p11/((PI4K3Ba + m31)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11));  // dwdp[17]
    dflux_v14_v_13_dp31 = 1.0*CERTERa*PI4K3Ba/(PI4K3Ba + m31);  // dwdp[18]
    dflux_v15_v_14_dp32 = 1.0*CERT;  // dwdp[19]
    dflux_v16_v_15_dp33 = 1.0*CERTTGNa*PKDDAGa/(PKDDAGa + m33);  // dwdp[20]
    dflux_v11_v_10_dpu3 = 1.0*u3;  // dwdp[21]
    dflux_v18_v_17_dpu4 = 1.0*u4;  // dwdp[22]
    dflux_v2_v_1_dpu5 = 1.0*PKD*p12*u5;  // dwdp[23]
    dflux_v3_v_2_dpu6 = 1.0*PKDDAGa*p13*u6;  // dwdp[24]
}

} // namespace model_Weber_BMC2015
} // namespace amici
