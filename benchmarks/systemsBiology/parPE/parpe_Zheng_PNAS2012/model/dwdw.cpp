#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 62> dwdw_colptrs_Zheng_PNAS2012_ = {
    0, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16
};

void dwdw_colptrs_Zheng_PNAS2012(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Zheng_PNAS2012_));
}
} // namespace model_Zheng_PNAS2012
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Zheng_PNAS2012 {

static constexpr std::array<sunindextype, 16> dwdw_rowvals_Zheng_PNAS2012_ = {
    1, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60
};

void dwdw_rowvals_Zheng_PNAS2012(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Zheng_PNAS2012_));
}
} // namespace model_Zheng_PNAS2012
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
namespace model_Zheng_PNAS2012 {

void dwdw_Zheng_PNAS2012(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_v1_v_0_dinflow = 1;  // dwdw[0]
    }

    // dynamic expressions
    dflux_v46_v_45_dinflow = K27me0K36me0;  // dwdw[1]
    dflux_v47_v_46_dinflow = K27me0K36me1;  // dwdw[2]
    dflux_v48_v_47_dinflow = K27me1K36me0;  // dwdw[3]
    dflux_v49_v_48_dinflow = K27me0K36me2;  // dwdw[4]
    dflux_v50_v_49_dinflow = K27me1K36me1;  // dwdw[5]
    dflux_v51_v_50_dinflow = K27me2K36me0;  // dwdw[6]
    dflux_v52_v_51_dinflow = K27me0K36me3;  // dwdw[7]
    dflux_v53_v_52_dinflow = K27me1K36me2;  // dwdw[8]
    dflux_v54_v_53_dinflow = K27me2K36me1;  // dwdw[9]
    dflux_v55_v_54_dinflow = K27me3K36me0;  // dwdw[10]
    dflux_v56_v_55_dinflow = K27me1K36me3;  // dwdw[11]
    dflux_v57_v_56_dinflow = K27me2K36me2;  // dwdw[12]
    dflux_v58_v_57_dinflow = K27me3K36me1;  // dwdw[13]
    dflux_v59_v_58_dinflow = K27me2K36me3;  // dwdw[14]
    dflux_v60_v_59_dinflow = K27me3K36me2;  // dwdw[15]
}

} // namespace model_Zheng_PNAS2012
} // namespace amici
