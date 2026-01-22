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

namespace amici {
namespace model_Weber_BMC2015 {

void w_Weber_BMC2015(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        flux_v4_v_3 = 1.0*s12;  // w[8]
        flux_v5_v_4 = 1.0*pu2*u2;  // w[9]
        flux_v10_v_9 = 1.0*s21;  // w[14]
        flux_v17_v_16 = 1.0*s31;  // w[21]
    }

    // dynamic expressions
    u3 = Ect_Expr_PI4K3beta_flag*Heaviside_0;  // w[1]
    u4 = Ect_Expr_CERT_flag*Heaviside_0;  // w[2]
    u5 = Heaviside_2*PdBu_dose;  // w[3]
    u6 = Heaviside_4*kb_NB142_70_dose;  // w[4]
    flux_v1_v_0 = 1.0*CERTERa*PI4K3Ba*PKD*p11*p31/((PI4K3Ba + m31)*(CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31) + m11));  // w[5]
    flux_v2_v_1 = 1.0*PKD*p12*(pu5*u5 + 1);  // w[6]
    flux_v3_v_2 = 1.0*PKDDAGa*p13*(pu6*u6 + 1);  // w[7]
    flux_v6_v_5 = 1.0*PKD*a11;  // w[10]
    flux_v7_v_6 = 1.0*PKDDAGa*a12;  // w[11]
    flux_v8_v_7 = 1.0*PI4K3Ba*p21;  // w[12]
    flux_v9_v_8 = 1.0*PI4K3B*PKDDAGa*p22/(PKDDAGa + m22);  // w[13]
    flux_v11_v_10 = 1.0*pu3*u3;  // w[15]
    flux_v12_v_11 = 1.0*PI4K3B*a21;  // w[16]
    flux_v13_v_12 = 1.0*PI4K3Ba*a22;  // w[17]
    flux_v14_v_13 = 1.0*CERTERa*PI4K3Ba*p31/(PI4K3Ba + m31);  // w[18]
    flux_v15_v_14 = 1.0*CERT*p32;  // w[19]
    flux_v16_v_15 = 1.0*CERTTGNa*PKDDAGa*p33/(PKDDAGa + m33);  // w[20]
    flux_v18_v_17 = 1.0*pu4*u4;  // w[22]
    flux_v19_v_18 = 1.0*CERTERa*a31;  // w[23]
    flux_v20_v_19 = 1.0*CERT*a32;  // w[24]
    flux_v21_v_20 = 1.0*CERTTGNa*a33;  // w[25]
}

} // namespace model_Weber_BMC2015
} // namespace amici
