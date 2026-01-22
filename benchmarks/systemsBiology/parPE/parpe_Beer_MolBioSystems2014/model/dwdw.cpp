#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 6> dwdw_colptrs_Beer_MolBioSystems2014_ = {
    0, 1, 1, 1, 1, 1
};

void dwdw_colptrs_Beer_MolBioSystems2014(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Beer_MolBioSystems2014_));
}
} // namespace model_Beer_MolBioSystems2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 1> dwdw_rowvals_Beer_MolBioSystems2014_ = {
    1
};

void dwdw_rowvals_Beer_MolBioSystems2014(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Beer_MolBioSystems2014_));
}
} // namespace model_Beer_MolBioSystems2014
} // namespace amici




#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "h.h"
#include "w.h"
#include "dwdw.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void dwdw_Beer_MolBioSystems2014(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_reaction1_dlag = 1.0*Bac*beta*(-Bac + Bacmax)/Bacmax;  // dwdw[0]
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
