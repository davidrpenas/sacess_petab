#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 5> dwdx_colptrs_Beer_MolBioSystems2014_ = {
    0, 3, 5, 7, 8
};

void dwdx_colptrs_Beer_MolBioSystems2014(SUNMatrixWrapper &dwdx){
    dwdx.set_indexptrs(gsl::make_span(dwdx_colptrs_Beer_MolBioSystems2014_));
}
} // namespace model_Beer_MolBioSystems2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 8> dwdx_rowvals_Beer_MolBioSystems2014_ = {
    1, 2, 3, 2, 3, 3, 4, 4
};

void dwdx_rowvals_Beer_MolBioSystems2014(SUNMatrixWrapper &dwdx){
    dwdx.set_indexvals(gsl::make_span(dwdx_rowvals_Beer_MolBioSystems2014_));
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
#include "dwdx.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void dwdx_Beer_MolBioSystems2014(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *spl, bool include_static){
    // static expressions
    if (include_static) {
        dflux_reaction4_dInd = -1.0*kdegi;  // dwdx[7]
    }

    // dynamic expressions
    dflux_reaction1_dBac = -1.0*Bac*beta*lag/Bacmax + 1.0*beta*lag*(-Bac + Bacmax)/Bacmax;  // dwdx[0]
    dflux_reaction2_dBac = -1.0*Glu*ksyn;  // dwdx[1]
    dflux_reaction3_dBac = 1.0*Glu*ksyn;  // dwdx[2]
    dflux_reaction2_dGlu = -1.0*Bac*ksyn;  // dwdx[3]
    dflux_reaction3_dGlu = 1.0*Bac*ksyn;  // dwdx[4]
    dflux_reaction3_dcGlu = -2.0*cGlu*kdim;  // dwdx[5]
    dflux_reaction4_dcGlu = 2.0*cGlu*kdim;  // dwdx[6]
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
