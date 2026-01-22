#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 13> dxdotdw_colptrs_Brannmark_JBC2010_ = {
    0, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22
};

void dxdotdw_colptrs_Brannmark_JBC2010(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Brannmark_JBC2010 {

static constexpr std::array<sunindextype, 22> dxdotdw_rowvals_Brannmark_JBC2010_ = {
    0, 1, 0, 1, 1, 2, 2, 3, 3, 4, 0, 2, 0, 4, 5, 6, 5, 6, 7, 8, 7, 8
};

void dxdotdw_rowvals_Brannmark_JBC2010(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Brannmark_JBC2010_));
}
} // namespace model_Brannmark_JBC2010
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
#include "dxdotdw.h"

namespace amici {
namespace model_Brannmark_JBC2010 {

void dxdotdw_Brannmark_JBC2010(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddIRdt_dflux_v1_v_0 = -1;  // dxdotdw[0]
    ddIRinsdt_dflux_v1_v_0 = 1;  // dxdotdw[1]
    ddIRdt_dflux_v2_v_1 = 1;  // dxdotdw[2]
    ddIRinsdt_dflux_v2_v_1 = -1;  // dxdotdw[3]
    ddIRinsdt_dflux_v3_v_2 = -1;  // dxdotdw[4]
    ddIRpdt_dflux_v3_v_2 = 1;  // dxdotdw[5]
    ddIRpdt_dflux_v4_v_3 = -1;  // dxdotdw[6]
    ddIRiPdt_dflux_v4_v_3 = 1;  // dxdotdw[7]
    ddIRiPdt_dflux_v5_v_4 = -1;  // dxdotdw[8]
    ddIRidt_dflux_v5_v_4 = 1;  // dxdotdw[9]
    ddIRdt_dflux_v6_v_5 = 1;  // dxdotdw[10]
    ddIRpdt_dflux_v6_v_5 = -1;  // dxdotdw[11]
    ddIRdt_dflux_v7_v_6 = 1;  // dxdotdw[12]
    ddIRidt_dflux_v7_v_6 = -1;  // dxdotdw[13]
    ddIRSdt_dflux_v8_v_7 = -1;  // dxdotdw[14]
    ddIRSiPdt_dflux_v8_v_7 = 1;  // dxdotdw[15]
    ddIRSdt_dflux_v9_v_8 = 1;  // dxdotdw[16]
    ddIRSiPdt_dflux_v9_v_8 = -1;  // dxdotdw[17]
    ddXdt_dflux_v10_v_9 = -1;  // dxdotdw[18]
    ddXpdt_dflux_v10_v_9 = 1;  // dxdotdw[19]
    ddXdt_dflux_v11_v_10 = 1;  // dxdotdw[20]
    ddXpdt_dflux_v11_v_10 = -1;  // dxdotdw[21]
}

} // namespace model_Brannmark_JBC2010
} // namespace amici
