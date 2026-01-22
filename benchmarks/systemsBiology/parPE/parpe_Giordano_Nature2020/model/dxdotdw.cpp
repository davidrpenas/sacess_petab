#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 31> dxdotdw_colptrs_Giordano_Nature2020_ = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 7, 9, 11, 14, 17, 19, 22, 24, 27, 29, 32
};

void dxdotdw_colptrs_Giordano_Nature2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Giordano_Nature2020_));
}
} // namespace model_Giordano_Nature2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Giordano_Nature2020 {

static constexpr std::array<sunindextype, 32> dxdotdw_rowvals_Giordano_Nature2020_ = {
    0, 1, 1, 2, 8, 1, 3, 1, 6, 2, 4, 2, 6, 9, 3, 4, 8, 3, 6, 3, 5, 8, 4, 5, 4, 6, 9, 5, 7, 5, 6, 9
};

void dxdotdw_rowvals_Giordano_Nature2020(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Giordano_Nature2020_));
}
} // namespace model_Giordano_Nature2020
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
namespace model_Giordano_Nature2020 {

void dxdotdw_Giordano_Nature2020(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSusceptibledt_dflux_Susceptible_to_Infected = -1;  // dxdotdw[0]
    ddInfecteddt_dflux_Susceptible_to_Infected = 1;  // dxdotdw[1]
    ddInfecteddt_dflux_Infected_to_Diagnosed = -1;  // dxdotdw[2]
    ddDiagnoseddt_dflux_Infected_to_Diagnosed = 1;  // dxdotdw[3]
    ddCumulativeDiagnoseddt_dflux_Infected_to_Diagnosed = 1;  // dxdotdw[4]
    ddInfecteddt_dflux_Infected_to_Ailing = -1;  // dxdotdw[5]
    ddAilingdt_dflux_Infected_to_Ailing = 1;  // dxdotdw[6]
    ddInfecteddt_dflux_Infected_to_Healed = -1;  // dxdotdw[7]
    ddHealeddt_dflux_Infected_to_Healed = 1;  // dxdotdw[8]
    ddDiagnoseddt_dflux_Diagnosed_to_Recognized = -1;  // dxdotdw[9]
    ddRecognizeddt_dflux_Diagnosed_to_Recognized = 1;  // dxdotdw[10]
    ddDiagnoseddt_dflux_Diagnosed_to_Healed = -1;  // dxdotdw[11]
    ddHealeddt_dflux_Diagnosed_to_Healed = 1;  // dxdotdw[12]
    ddDiagnosedHealeddt_dflux_Diagnosed_to_Healed = 1;  // dxdotdw[13]
    ddAilingdt_dflux_Ailing_to_Recognised = -1;  // dxdotdw[14]
    ddRecognizeddt_dflux_Ailing_to_Recognised = 1;  // dxdotdw[15]
    ddCumulativeDiagnoseddt_dflux_Ailing_to_Recognised = 1;  // dxdotdw[16]
    ddAilingdt_dflux_Ailing_to_Healed = -1;  // dxdotdw[17]
    ddHealeddt_dflux_Ailing_to_Healed = 1;  // dxdotdw[18]
    ddAilingdt_dflux_Ailing_to_Threatened = -1;  // dxdotdw[19]
    ddThreateneddt_dflux_Ailing_to_Threatened = 1;  // dxdotdw[20]
    ddCumulativeDiagnoseddt_dflux_Ailing_to_Threatened = 1;  // dxdotdw[21]
    ddRecognizeddt_dflux_Recognised_to_Threatened = -1;  // dxdotdw[22]
    ddThreateneddt_dflux_Recognised_to_Threatened = 1;  // dxdotdw[23]
    ddRecognizeddt_dflux_Recognised_to_Healed = -1;  // dxdotdw[24]
    ddHealeddt_dflux_Recognised_to_Healed = 1;  // dxdotdw[25]
    ddDiagnosedHealeddt_dflux_Recognised_to_Healed = 1;  // dxdotdw[26]
    ddThreateneddt_dflux_Threatened_to_Extinct = -1;  // dxdotdw[27]
    ddExtinctdt_dflux_Threatened_to_Extinct = 1;  // dxdotdw[28]
    ddThreateneddt_dflux_Threatened_to_Healed = -1;  // dxdotdw[29]
    ddHealeddt_dflux_Threatened_to_Healed = 1;  // dxdotdw[30]
    ddDiagnosedHealeddt_dflux_Threatened_to_Healed = 1;  // dxdotdw[31]
}

} // namespace model_Giordano_Nature2020
} // namespace amici
