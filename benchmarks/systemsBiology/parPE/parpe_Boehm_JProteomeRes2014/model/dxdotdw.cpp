#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<sunindextype, 11> dxdotdw_colptrs_Boehm_JProteomeRes2014_ = {
    0, 0, 2, 5, 7, 9, 11, 13, 15, 18, 20
};

void dxdotdw_colptrs_Boehm_JProteomeRes2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexptrs(gsl::make_span(dxdotdw_colptrs_Boehm_JProteomeRes2014_));
}
} // namespace model_Boehm_JProteomeRes2014
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

static constexpr std::array<sunindextype, 20> dxdotdw_rowvals_Boehm_JProteomeRes2014_ = {
    0, 3, 0, 1, 2, 1, 4, 3, 5, 2, 6, 4, 7, 0, 5, 0, 1, 6, 1, 7
};

void dxdotdw_rowvals_Boehm_JProteomeRes2014(SUNMatrixWrapper &dxdotdw){
    dxdotdw.set_indexvals(gsl::make_span(dxdotdw_rowvals_Boehm_JProteomeRes2014_));
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
#include "dxdotdw.h"

namespace amici {
namespace model_Boehm_JProteomeRes2014 {

void dxdotdw_Boehm_JProteomeRes2014(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    ddSTAT5Adt_dflux_v1_v_0 = -1.4285714285714286;  // dxdotdw[0]
    ddpApAdt_dflux_v1_v_0 = 0.7142857142857143;  // dxdotdw[1]
    ddSTAT5Adt_dflux_v2_v_1 = -0.7142857142857143;  // dxdotdw[2]
    ddSTAT5Bdt_dflux_v2_v_1 = -0.7142857142857143;  // dxdotdw[3]
    ddpApBdt_dflux_v2_v_1 = 0.7142857142857143;  // dxdotdw[4]
    ddSTAT5Bdt_dflux_v3_v_2 = -1.4285714285714286;  // dxdotdw[5]
    ddpBpBdt_dflux_v3_v_2 = 0.7142857142857143;  // dxdotdw[6]
    ddpApAdt_dflux_v4_v_3 = -0.7142857142857143;  // dxdotdw[7]
    ddnucpApAdt_dflux_v4_v_3 = 2.2222222222222223;  // dxdotdw[8]
    ddpApBdt_dflux_v5_v_4 = -0.7142857142857143;  // dxdotdw[9]
    ddnucpApBdt_dflux_v5_v_4 = 2.2222222222222223;  // dxdotdw[10]
    ddpBpBdt_dflux_v6_v_5 = -0.7142857142857143;  // dxdotdw[11]
    ddnucpBpBdt_dflux_v6_v_5 = 2.2222222222222223;  // dxdotdw[12]
    ddSTAT5Adt_dflux_v7_v_6 = 1.4285714285714286;  // dxdotdw[13]
    ddnucpApAdt_dflux_v7_v_6 = -2.2222222222222223;  // dxdotdw[14]
    ddSTAT5Adt_dflux_v8_v_7 = 0.7142857142857143;  // dxdotdw[15]
    ddSTAT5Bdt_dflux_v8_v_7 = 0.7142857142857143;  // dxdotdw[16]
    ddnucpApBdt_dflux_v8_v_7 = -2.2222222222222223;  // dxdotdw[17]
    ddSTAT5Bdt_dflux_v9_v_8 = 1.4285714285714286;  // dxdotdw[18]
    ddnucpBpBdt_dflux_v9_v_8 = -2.2222222222222223;  // dxdotdw[19]
}

} // namespace model_Boehm_JProteomeRes2014
} // namespace amici
