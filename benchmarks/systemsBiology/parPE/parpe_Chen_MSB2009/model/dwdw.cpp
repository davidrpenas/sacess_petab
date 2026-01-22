#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Chen_MSB2009 {

static constexpr std::array<sunindextype, 836> dwdw_colptrs_Chen_MSB2009_ = {
    0, 30, 36, 40, 55, 56, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86
};

void dwdw_colptrs_Chen_MSB2009(SUNMatrixWrapper &dwdw){
    dwdw.set_indexptrs(gsl::make_span(dwdw_colptrs_Chen_MSB2009_));
}
} // namespace model_Chen_MSB2009
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_Chen_MSB2009 {

static constexpr std::array<sunindextype, 86> dwdw_rowvals_Chen_MSB2009_ = {
    28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 794, 797, 820, 821, 827, 828, 829, 830, 8, 11, 777, 778, 788, 791, 8, 11, 788, 791, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 179, 792, 762, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 794, 797, 820, 821, 827, 828, 829, 830
};

void dwdw_rowvals_Chen_MSB2009(SUNMatrixWrapper &dwdw){
    dwdw.set_indexvals(gsl::make_span(dwdw_rowvals_Chen_MSB2009_));
}
} // namespace model_Chen_MSB2009
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
namespace model_Chen_MSB2009 {

void dwdw_Chen_MSB2009(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, bool include_static){

    // dynamic expressions
    dflux_v21_v21_dStep1 = -1.0*c123;  // dwdw[0]
    dflux_v22_v22_dStep1 = -1.0*c124;  // dwdw[1]
    dflux_v23_v23_dStep1 = -1.0*c125;  // dwdw[2]
    dflux_v24_v24_dStep1 = -1.0*c126;  // dwdw[3]
    dflux_v25_v25_dStep1 = -1.0*c116;  // dwdw[4]
    dflux_v26_v26_dStep1 = -1.0*c122;  // dwdw[5]
    dflux_v27_v27_dStep1 = -1.0*c127;  // dwdw[6]
    dflux_v28_v28_dStep1 = -1.0*c128;  // dwdw[7]
    dflux_v29_v29_dStep1 = -1.0*c129;  // dwdw[8]
    dflux_v30_v30_dStep1 = -1.0*c130;  // dwdw[9]
    dflux_v31_v31_dStep1 = -1.0*c131;  // dwdw[10]
    dflux_v32_v32_dStep1 = -1.0*c132;  // dwdw[11]
    dflux_v33_v33_dStep1 = -1.0*c133;  // dwdw[12]
    dflux_v34_v34_dStep1 = -1.0*c134;  // dwdw[13]
    dflux_v35_v35_dStep1 = -1.0*c135;  // dwdw[14]
    dflux_v36_v36_dStep1 = -1.0*c136;  // dwdw[15]
    dflux_v37_v37_dStep1 = -1.0*c137;  // dwdw[16]
    dflux_v38_v38_dStep1 = -1.0*c138;  // dwdw[17]
    dflux_v39_v39_dStep1 = -1.0*c139;  // dwdw[18]
    dflux_v40_v40_dStep1 = -1.0*c168;  // dwdw[19]
    dflux_v41_v41_dStep1 = -1.0*c169;  // dwdw[20]
    dflux_v42_v42_dStep1 = -1.0*c170;  // dwdw[21]
    dflux_v787_v802_dStep1 = -1.0*c115;  // dwdw[22]
    dflux_v790_v805_dStep1 = -1.0*c121;  // dwdw[23]
    dflux_v813_v828_dStep1 = -1.0*c2;  // dwdw[24]
    dflux_v814_v829_dStep1 = -1.0*c524;  // dwdw[25]
    dflux_v820_v855_dStep1 = -1.0*c555;  // dwdw[26]
    dflux_v821_v856_dStep1 = -1.0*c556;  // dwdw[27]
    dflux_v822_v857_dStep1 = -1.0*c557;  // dwdw[28]
    dflux_v823_v858_dStep1 = -1.0*c558;  // dwdw[29]
    dflux_v1_v1_dStep10Mio = 1.0*c1*c2;  // dwdw[30]
    dflux_v4_v4_dStep10Mio = 1.0*c1*c286;  // dwdw[31]
    dflux_v770_v784_dStep10Mio = 1.0*c140*c514;  // dwdw[32]
    dflux_v771_v785_dStep10Mio = 1.0*c143*c514;  // dwdw[33]
    dflux_v781_v795_dStep10Mio = 1.0*c1*c525;  // dwdw[34]
    dflux_v784_v798_dStep10Mio = 1.0*c1*c524;  // dwdw[35]
    dflux_v1_v1_dStep33 = -1.0*c3;  // dwdw[36]
    dflux_v4_v4_dStep33 = -1.0*c499;  // dwdw[37]
    dflux_v781_v795_dStep33 = -1.0*c526;  // dwdw[38]
    dflux_v784_v798_dStep33 = -1.0*c529;  // dwdw[39]
    dflux_v155_v163_dStep5 = -1.0*c20;  // dwdw[40]
    dflux_v156_v164_dStep5 = -1.0*c6;  // dwdw[41]
    dflux_v157_v165_dStep5 = -1.0*c8;  // dwdw[42]
    dflux_v158_v166_dStep5 = -1.0*c17;  // dwdw[43]
    dflux_v159_v167_dStep5 = -1.0*c63;  // dwdw[44]
    dflux_v160_v168_dStep5 = -1.0*c64;  // dwdw[45]
    dflux_v161_v169_dStep5 = -1.0*c19;  // dwdw[46]
    dflux_v162_v170_dStep5 = -1.0*c21;  // dwdw[47]
    dflux_v163_v171_dStep5 = -1.0*c65;  // dwdw[48]
    dflux_v164_v172_dStep5 = -1.0*c66;  // dwdw[49]
    dflux_v165_v173_dStep5 = -1.0*c67;  // dwdw[50]
    dflux_v166_v174_dStep5 = -1.0*c68;  // dwdw[51]
    dflux_v167_v175_dStep5 = -1.0*c18;  // dwdw[52]
    dflux_v172_v180_dStep5 = -1.0*c350;  // dwdw[53]
    dflux_v785_v799_dStep5 = -1.0*c530;  // dwdw[54]
    dflux_v755_v769_dStepLate = 1.0*c60;  // dwdw[55]
    dflux_v21_v21_dStepMini = 1.0*c105*c159;  // dwdw[56]
    dflux_v22_v22_dStepMini = 1.0*c105*c160;  // dwdw[57]
    dflux_v23_v23_dStepMini = 1.0*c105*c161;  // dwdw[58]
    dflux_v24_v24_dStepMini = 1.0*c105*c11;  // dwdw[59]
    dflux_v25_v25_dStepMini = 1.0*c105*c4;  // dwdw[60]
    dflux_v26_v26_dStepMini = 1.0*c105*c145;  // dwdw[61]
    dflux_v27_v27_dStepMini = 1.0*c105*c146;  // dwdw[62]
    dflux_v28_v28_dStepMini = 1.0*c105*c147;  // dwdw[63]
    dflux_v29_v29_dStepMini = 1.0*c105*c284;  // dwdw[64]
    dflux_v30_v30_dStepMini = 1.0*c105*c427;  // dwdw[65]
    dflux_v31_v31_dStepMini = 1.0*c105*c428;  // dwdw[66]
    dflux_v32_v32_dStepMini = 1.0*c105*c429;  // dwdw[67]
    dflux_v33_v33_dStepMini = 1.0*c105*c436;  // dwdw[68]
    dflux_v34_v34_dStepMini = 1.0*c105*c439;  // dwdw[69]
    dflux_v35_v35_dStepMini = 1.0*c105*c442;  // dwdw[70]
    dflux_v36_v36_dStepMini = 1.0*c105*c483;  // dwdw[71]
    dflux_v37_v37_dStepMini = 1.0*c105*c516;  // dwdw[72]
    dflux_v38_v38_dStepMini = 1.0*c105*c517;  // dwdw[73]
    dflux_v39_v39_dStepMini = 1.0*c105*c345;  // dwdw[74]
    dflux_v40_v40_dStepMini = 1.0*c105*c355;  // dwdw[75]
    dflux_v41_v41_dStepMini = 1.0*c105*c421;  // dwdw[76]
    dflux_v42_v42_dStepMini = 1.0*c105*c422;  // dwdw[77]
    dflux_v787_v802_dStepMini = 1.0*c105*c500;  // dwdw[78]
    dflux_v790_v805_dStepMini = 1.0*c105*c527;  // dwdw[79]
    dflux_v813_v828_dStepMini = 1.0*c105*c531;  // dwdw[80]
    dflux_v814_v829_dStepMini = 1.0*c105*c532;  // dwdw[81]
    dflux_v820_v855_dStepMini = 1.0*c105*c550;  // dwdw[82]
    dflux_v821_v856_dStepMini = 1.0*c105*c551;  // dwdw[83]
    dflux_v822_v857_dStepMini = 1.0*c105*c552;  // dwdw[84]
    dflux_v823_v858_dStepMini = 1.0*c105*c553;  // dwdw[85]
}

} // namespace model_Chen_MSB2009
} // namespace amici
