#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 27> dwdw_colptrs_Weber_BMC2015_ = {
    0, 1, 2, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
};

void dwdw_colptrs_Weber_BMC2015(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Weber_BMC2015_));
}
} // namespace model_Weber_BMC2015
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Weber_BMC2015 {

static constexpr std::array<sunindextype, 5> dwdw_rowvals_Weber_BMC2015_ = {
    9, 15, 22, 6, 7
};

void dwdw_rowvals_Weber_BMC2015(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Weber_BMC2015_));
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
#include "dwdw.h"

namespace amici {
namespace model_Weber_BMC2015 {

void dwdw_Weber_BMC2015(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v5_v_4_du2 = 1.0*pu2;  // dwdw[0]
        dflux_v11_v_10_du3 = 1.0*pu3;  // dwdw[1]
        dflux_v18_v_17_du4 = 1.0*pu4;  // dwdw[2]
    }

    // dynamic expressions
    dflux_v2_v_1_du5 = 1.0*PKD*p12*pu5;  // dwdw[3]
    dflux_v3_v_2_du6 = 1.0*PKDDAGa*p13*pu6;  // dwdw[4]
}

} // namespace model_Weber_BMC2015
} // namespace amici
