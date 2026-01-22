#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 31> dwdp_colptrs_Laske_PLOSComputBiol2019_ = {
    0, 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6
};

void dwdp_colptrs_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Laske_PLOSComputBiol2019_));
}
} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Laske_PLOSComputBiol2019 {

static constexpr std::array<sunindextype, 6> dwdp_rowvals_Laske_PLOSComputBiol2019_ = {
    61, 62, 98, 16, 100, 67
};

void dwdp_rowvals_Laske_PLOSComputBiol2019(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Laske_PLOSComputBiol2019_));
}
} // namespace model_Laske_PLOSComputBiol2019
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
namespace model_Laske_PLOSComputBiol2019 {

void dwdp_Laske_PLOSComputBiol2019(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){
    // static expressions
    if (include_static) {
        dModelValue_79_dk_syn_R_M = 1;  // dwdp[3]
    }

    // dynamic expressions
    dflux_V_cyt_Import_dk_imp = 1.0*Vp_cyt;  // dwdp[0]
    dflux_R_C_Synthesis_dk_syn_R_C = 1.0*Vp_nuc;  // dwdp[1]
    dflux_Vp_nuc_binding_of_P_M1_dk_bind_M1 = 1.0*P_M1*Vp_nuc;  // dwdp[2]
    dflux_Vp_cyt_M1_Virion_release_dk_rel = 1.0*P_HA*P_M1*P_M2*P_NA*P_NEP*P_NP*P_RdRp*Vp_cyt_M1/((KmB + P_RdRp)*(KmC + P_HA)*(KmD + P_NP)*(KmE + P_NA)*(KmF + P_M1)*(KmG + P_M2)*(KmH + P_NEP));  // dwdp[4]
    dflux_R_V_synthesis_dk_syn_R_V = 1.0*Cp;  // dwdp[5]
}

} // namespace model_Laske_PLOSComputBiol2019
} // namespace amici
