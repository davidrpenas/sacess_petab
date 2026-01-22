#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<sunindextype, 9> dxdotdw_colptrs_Armistead_CellDeathDis2024_ = {
    0, 1, 3, 5, 7, 8, 10, 12, 13
};

void dxdotdw_colptrs_Armistead_CellDeathDis2024(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Armistead_CellDeathDis2024_));
}
} // namespace model_Armistead_CellDeathDis2024
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

static constexpr std::array<sunindextype, 13> dxdotdw_rowvals_Armistead_CellDeathDis2024_ = {
    3, 0, 3, 0, 1, 0, 1, 0, 1, 2, 1, 2, 2
};

void dxdotdw_rowvals_Armistead_CellDeathDis2024(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Armistead_CellDeathDis2024_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Armistead_CellDeathDis2024 {

void dxdotdw_Armistead_CellDeathDis2024(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSphingadt_dflux_R1 = 1;  // dxdotdw[0]
    ddCerdt_dflux_R2 = 1;  // dxdotdw[1]
    ddSphingadt_dflux_R2 = -1;  // dxdotdw[2]
    ddCerdt_dflux_R3 = -1;  // dxdotdw[3]
    ddSphingodt_dflux_R3 = 1;  // dxdotdw[4]
    ddCerdt_dflux_R4 = 1;  // dxdotdw[5]
    ddSphingodt_dflux_R4 = -1;  // dxdotdw[6]
    ddCerdt_dflux_R5 = -1;  // dxdotdw[7]
    ddSphingodt_dflux_R6 = -1;  // dxdotdw[8]
    ddS1Pdt_dflux_R6 = 1;  // dxdotdw[9]
    ddSphingodt_dflux_R7 = 1;  // dxdotdw[10]
    ddS1Pdt_dflux_R7 = -1;  // dxdotdw[11]
    ddS1Pdt_dflux_R8 = -1;  // dxdotdw[12]
}

} // namespace model_Armistead_CellDeathDis2024
} // namespace amici
