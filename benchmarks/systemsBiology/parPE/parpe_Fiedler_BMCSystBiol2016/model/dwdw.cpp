#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 8> dwdw_colptrs_Fiedler_BMCSystBiol2016_ = {
    0, 1, 1, 1, 1, 1, 1, 1
};

void dwdw_colptrs_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Fiedler_BMCSystBiol2016_));
}
} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

static constexpr std::array<sunindextype, 1> dwdw_rowvals_Fiedler_BMCSystBiol2016_ = {
    1
};

void dwdw_rowvals_Fiedler_BMCSystBiol2016(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Fiedler_BMCSystBiol2016_));
}
} // namespace model_Fiedler_BMCSystBiol2016
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
#include "dwdw.h"

namespace amici {
namespace model_Fiedler_BMCSystBiol2016 {

void dwdw_Fiedler_BMCSystBiol2016(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_v1_v_0_dk1max = 1.0*K_1*RAF/(K_1 + pERK);  // dwdw[0]
}

} // namespace model_Fiedler_BMCSystBiol2016
} // namespace amici
