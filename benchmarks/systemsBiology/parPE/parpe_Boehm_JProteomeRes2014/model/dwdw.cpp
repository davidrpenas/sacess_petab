#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<sunindextype, 11> dwdw_colptrs_Boehm_JProteomeRes2014_ = {
    0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
};

void dwdw_colptrs_Boehm_JProteomeRes2014(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Boehm_JProteomeRes2014_));
}
} // namespace model_Boehm_JProteomeRes2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<sunindextype, 3> dwdw_rowvals_Boehm_JProteomeRes2014_ = {
    1, 2, 3
};

void dwdw_rowvals_Boehm_JProteomeRes2014(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Boehm_JProteomeRes2014_));
}
} // namespace model_Boehm_JProteomeRes2014
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
namespace model_Boehm_JProteomeRes2014 {

void dwdw_Boehm_JProteomeRes2014(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_v1_v_0_dBaF3_Epo = 1.3999999999999999*std::pow(STAT5A, 2)*k_phos;  // dwdw[0]
    dflux_v2_v_1_dBaF3_Epo = 1.3999999999999999*STAT5A*STAT5B*k_phos;  // dwdw[1]
    dflux_v3_v_2_dBaF3_Epo = 1.3999999999999999*std::pow(STAT5B, 2)*k_phos;  // dwdw[2]
}

} // namespace model_Boehm_JProteomeRes2014
} // namespace amici
