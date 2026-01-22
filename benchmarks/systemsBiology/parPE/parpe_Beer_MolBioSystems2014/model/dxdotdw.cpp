#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 6> dxdotdw_colptrs_Beer_MolBioSystems2014_ = {
    0, 0, 1, 2, 3, 4
};

void dxdotdw_colptrs_Beer_MolBioSystems2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Beer_MolBioSystems2014_));
}
} // namespace model_Beer_MolBioSystems2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Beer_MolBioSystems2014 {

static constexpr std::array<sunindextype, 4> dxdotdw_rowvals_Beer_MolBioSystems2014_ = {
    0, 1, 2, 3
};

void dxdotdw_rowvals_Beer_MolBioSystems2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Beer_MolBioSystems2014_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Beer_MolBioSystems2014 {

void dxdotdw_Beer_MolBioSystems2014(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddBacdt_dflux_reaction1 = 1;  // dxdotdw[0]
    ddGludt_dflux_reaction2 = 1;  // dxdotdw[1]
    ddcGludt_dflux_reaction3 = 1;  // dxdotdw[2]
    ddInddt_dflux_reaction4 = 1;  // dxdotdw[3]
}

} // namespace model_Beer_MolBioSystems2014
} // namespace amici
