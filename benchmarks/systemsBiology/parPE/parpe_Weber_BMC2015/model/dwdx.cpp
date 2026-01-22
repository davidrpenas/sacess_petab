#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 8> dwdx_colptrs_Weber_BMC2015_ = {
    0, 3, 7, 9, 13, 16, 18, 20
};

void dwdx_colptrs_Weber_BMC2015(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Weber_BMC2015_));
}
} // namespace model_Weber_BMC2015
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 20> dwdx_rowvals_Weber_BMC2015_ = {
    5, 6, 10, 7, 11, 13, 20, 13, 16, 5, 12, 17, 18, 5, 18, 23, 19, 24, 20, 25
};

void dwdx_rowvals_Weber_BMC2015(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Weber_BMC2015_));
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
#include "dwdx.h"

namespace amici {
namespace model_Weber_BMC2015 {

void dwdx_Weber_BMC2015(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v6_v_5_dPKD = 1.0*a11;  // dwdx[2]
        dflux_v7_v_6_dPKDDAGa = 1.0*a12;  // dwdx[4]
        dflux_v12_v_11_dPI4K3B = 1.0*a21;  // dwdx[8]
        dflux_v8_v_7_dPI4K3Ba = 1.0*p21;  // dwdx[10]
        dflux_v13_v_12_dPI4K3Ba = 1.0*a22;  // dwdx[11]
        dflux_v19_v_18_dCERTERa = 1.0*a31;  // dwdx[15]
        dflux_v15_v_14_dCERT = 1.0*p32;  // dwdx[16]
        dflux_v20_v_19_dCERT = 1.0*a32;  // dwdx[17]
        dflux_v21_v_20_dCERTTGNa = 1.0*a33;  // dwdx[19]
    }

    // dynamic expressions
    dflux_v1_v_0_dPKD = 1.0*CERTERa*PI4K3Ba*p11*p31/((PI4K3Ba + m31)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11));  // dwdx[0]
    dflux_v2_v_1_dPKD = 1.0*p12*(pu5*u5 + 1);  // dwdx[1]
    dflux_v3_v_2_dPKDDAGa = 1.0*p13*(pu6*u6 + 1);  // dwdx[3]
    dflux_v9_v_8_dPKDDAGa = -1.0*PI4K3B*PKDDAGa*p22/std::pow(PKDDAGa + m22, 2) + 1.0*PI4K3B*p22/(PKDDAGa + m22);  // dwdx[5]
    dflux_v16_v_15_dPKDDAGa = -1.0*CERTTGNa*PKDDAGa*p33/std::pow(PKDDAGa + m33, 2) + 1.0*CERTTGNa*p33/(PKDDAGa + m33);  // dwdx[6]
    dflux_v9_v_8_dPI4K3B = 1.0*PKDDAGa*p22/(PKDDAGa + m22);  // dwdx[7]
    dflux_v1_v_0_dPI4K3Ba = -1.0*CERTERa*PI4K3Ba*PKD*p11*p31*(-CERTERa*PI4K3Ba*p31/std::pow(PI4K3Ba + m31, 2) + CERTERa*p31/(PI4K3Ba + m31))/((PI4K3Ba + m31)*std::pow(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11, 2)) - 1.0*CERTERa*PI4K3Ba*PKD*p11*p31/(std::pow(PI4K3Ba + m31, 2)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11)) + 1.0*CERTERa*PKD*p11*p31/((PI4K3Ba + m31)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11));  // dwdx[9]
    dflux_v14_v_13_dPI4K3Ba = -1.0*CERTERa*PI4K3Ba*p31/std::pow(PI4K3Ba + m31, 2) + 1.0*CERTERa*p31/(PI4K3Ba + m31);  // dwdx[12]
    dflux_v1_v_0_dCERTERa = -1.0*CERTERa*std::pow(PI4K3Ba, 2)*PKD*p11*std::pow(p31, 2)/(std::pow(PI4K3Ba + m31, 2)*std::pow(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11, 2)) + 1.0*PI4K3Ba*PKD*p11*p31/((PI4K3Ba + m31)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11));  // dwdx[13]
    dflux_v14_v_13_dCERTERa = 1.0*PI4K3Ba*p31/(PI4K3Ba + m31);  // dwdx[14]
    dflux_v16_v_15_dCERTTGNa = 1.0*PKDDAGa*p33/(PKDDAGa + m33);  // dwdx[18]
}

} // namespace model_Weber_BMC2015
} // namespace amici
