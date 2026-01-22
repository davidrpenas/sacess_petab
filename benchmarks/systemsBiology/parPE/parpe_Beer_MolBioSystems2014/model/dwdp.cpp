#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 10> dwdp_colptrs_Beer_MolBioSystems2014_ = {
    0, 1, 2, 3, 5, 7, 7, 7, 7, 7
};

void dwdp_colptrs_Beer_MolBioSystems2014(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_Beer_MolBioSystems2014_));
}
} // namespace model_Beer_MolBioSystems2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 7> dwdp_rowvals_Beer_MolBioSystems2014_ = {
    1, 1, 4, 3, 4, 2, 3
};

void dwdp_rowvals_Beer_MolBioSystems2014(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_Beer_MolBioSystems2014_));
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
#include "dwdp.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void dwdp_Beer_MolBioSystems2014(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_reaction1_dBacmax = 1.0*Bac*beta*lag/Bacmax - 1.0*Bac*beta*lag*(-Bac + Bacmax)/std::pow(Bacmax, 2);  // dwdp[0]
    dflux_reaction1_dbeta = 1.0*Bac*lag*(-Bac + Bacmax)/Bacmax;  // dwdp[1]
    dflux_reaction4_dkdegi = -1.0*Ind;  // dwdp[2]
    dflux_reaction3_dkdim = -1.0*std::pow(cGlu, 2);  // dwdp[3]
    dflux_reaction4_dkdim = 1.0*std::pow(cGlu, 2);  // dwdp[4]
    dflux_reaction2_dksyn = -1.0*Bac*Glu;  // dwdp[5]
    dflux_reaction3_dksyn = 1.0*Bac*Glu;  // dwdp[6]
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
