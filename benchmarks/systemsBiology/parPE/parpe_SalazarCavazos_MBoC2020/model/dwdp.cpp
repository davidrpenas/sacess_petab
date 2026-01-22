#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

static constexpr std::array<sunindextype, 25> dwdp_colptrs_SalazarCavazos_MBoC2020_ = {
    0, 0, 0, 72, 126, 150, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168
};

void dwdp_colptrs_SalazarCavazos_MBoC2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexptrs(gsl::make_span(dwdp_colptrs_SalazarCavazos_MBoC2020_));
}
} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici

#include "amici/sundials_matrix_wrapper.h"
#include "sundials/sundials_types.h"

#include <array>
#include <algorithm>

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

static constexpr std::array<sunindextype, 168> dwdp_rowvals_SalazarCavazos_MBoC2020_ = {
    12, 13, 15, 16, 24, 25, 26, 27, 30, 31, 32, 33, 36, 37, 38, 39, 74, 75, 76, 77, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 200, 201, 204, 205, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 538, 539, 540, 541, 542, 543, 608, 609, 14, 17, 28, 29, 34, 35, 40, 41, 78, 79, 80, 85, 86, 87, 98, 99, 100, 101, 102, 202, 203, 206, 207, 222, 223, 224, 225, 226, 227, 228, 229, 230, 386, 387, 400, 401, 402, 403, 404, 405, 406, 407, 408, 544, 545, 546, 547, 548, 549, 550, 610, 611, 612, 622, 12, 13, 15, 16, 24, 25, 26, 27, 30, 31, 32, 33, 74, 75, 76, 77, 81, 82, 83, 84, 200, 201, 204, 205, 14, 17, 28, 29, 34, 35, 78, 79, 80, 85, 86, 87, 202, 203, 206, 207, 386, 387
};

void dwdp_rowvals_SalazarCavazos_MBoC2020(SUNMatrixWrapper &dwdp){
    dwdp.set_indexvals(gsl::make_span(dwdp_rowvals_SalazarCavazos_MBoC2020_));
}
} // namespace model_SalazarCavazos_MBoC2020
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
#include "dwdp.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void dwdp_SalazarCavazos_MBoC2020(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl, bool include_static){

    // dynamic expressions
    dflux_reaction_6_dkdephosY1068__FREE = kphos_f*ratio_kpkd_Y1068__FREE*species_6;  // dwdp[0]
    dflux_reaction_7_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kpkd_Y1068__FREE*species_6;  // dwdp[1]
    dflux_reaction_9_dkdephosY1068__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_7;  // dwdp[2]
    dflux_reaction_10_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_7;  // dwdp[3]
    dflux_reaction_18_dkdephosY1068__FREE = kphos_f*ratio_kpkd_Y1068__FREE*species_9;  // dwdp[4]
    dflux_reaction_19_dkdephosY1068__FREE = kphos_f*ratio_kpkd_Y1068__FREE*species_10;  // dwdp[5]
    dflux_reaction_20_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kpkd_Y1068__FREE*species_8;  // dwdp[6]
    dflux_reaction_21_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kpkd_Y1068__FREE*species_10;  // dwdp[7]
    dflux_reaction_24_dkdephosY1068__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_12;  // dwdp[8]
    dflux_reaction_25_dkdephosY1068__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_13;  // dwdp[9]
    dflux_reaction_26_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_11;  // dwdp[10]
    dflux_reaction_27_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_13;  // dwdp[11]
    dflux_reaction_30_dkdephosY1068__FREE = kdephosY1068_f*species_8;  // dwdp[12]
    dflux_reaction_31_dkdephosY1068__FREE = kdephosY1068_f*species_11;  // dwdp[13]
    dflux_reaction_32_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_9;  // dwdp[14]
    dflux_reaction_33_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_12;  // dwdp[15]
    dflux_reaction_68_dkdephosY1068__FREE = kphos_f*ratio_kpkd_Y1068__FREE*species_19;  // dwdp[16]
    dflux_reaction_69_dkdephosY1068__FREE = kphos_f*ratio_kpkd_Y1068__FREE*species_25;  // dwdp[17]
    dflux_reaction_70_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kpkd_Y1068__FREE*species_18;  // dwdp[18]
    dflux_reaction_71_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kpkd_Y1068__FREE*species_23;  // dwdp[19]
    dflux_reaction_75_dkdephosY1068__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_22;  // dwdp[20]
    dflux_reaction_76_dkdephosY1068__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_26;  // dwdp[21]
    dflux_reaction_77_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_21;  // dwdp[22]
    dflux_reaction_78_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_24;  // dwdp[23]
    dflux_reaction_82_dkdephosY1068__FREE = kdephosY1068_f*species_14;  // dwdp[24]
    dflux_reaction_83_dkdephosY1068__FREE = kdephosY1068_f*species_17;  // dwdp[25]
    dflux_reaction_84_dkdephosY1068__FREE = kdephosY1068_f*species_18;  // dwdp[26]
    dflux_reaction_85_dkdephosY1068__FREE = kdephosY1068_f*species_20;  // dwdp[27]
    dflux_reaction_86_dkdephosY1068__FREE = kdephosY1068_f*species_21;  // dwdp[28]
    dflux_reaction_87_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_15;  // dwdp[29]
    dflux_reaction_88_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_17;  // dwdp[30]
    dflux_reaction_89_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_19;  // dwdp[31]
    dflux_reaction_90_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_20;  // dwdp[32]
    dflux_reaction_91_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_22;  // dwdp[33]
    dflux_reaction_194_dkdephosY1068__FREE = kphos_f*ratio_kpkd_Y1068__FREE*species_39;  // dwdp[34]
    dflux_reaction_195_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kpkd_Y1068__FREE*species_38;  // dwdp[35]
    dflux_reaction_198_dkdephosY1068__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_44;  // dwdp[36]
    dflux_reaction_199_dkdephosY1068__FREE = kphos_f*ratio_kphosY1173*ratio_kphos_receiver*ratio_kpkd_Y1068__FREE*species_43;  // dwdp[37]
    dflux_reaction_202_dkdephosY1068__FREE = kdephosY1068_f*species_27;  // dwdp[38]
    dflux_reaction_203_dkdephosY1068__FREE = kdephosY1068_f*species_30;  // dwdp[39]
    dflux_reaction_204_dkdephosY1068__FREE = kdephosY1068_f*species_31;  // dwdp[40]
    dflux_reaction_205_dkdephosY1068__FREE = kdephosY1068_f*species_35;  // dwdp[41]
    dflux_reaction_206_dkdephosY1068__FREE = kdephosY1068_f*species_36;  // dwdp[42]
    dflux_reaction_207_dkdephosY1068__FREE = kdephosY1068_f*species_40;  // dwdp[43]
    dflux_reaction_208_dkdephosY1068__FREE = kdephosY1068_f*species_41;  // dwdp[44]
    dflux_reaction_209_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_28;  // dwdp[45]
    dflux_reaction_210_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_30;  // dwdp[46]
    dflux_reaction_211_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_32;  // dwdp[47]
    dflux_reaction_212_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_35;  // dwdp[48]
    dflux_reaction_213_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_37;  // dwdp[49]
    dflux_reaction_214_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_40;  // dwdp[50]
    dflux_reaction_215_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_42;  // dwdp[51]
    dflux_reaction_382_dkdephosY1068__FREE = kdephosY1068_f*species_45;  // dwdp[52]
    dflux_reaction_383_dkdephosY1068__FREE = kdephosY1068_f*species_46;  // dwdp[53]
    dflux_reaction_384_dkdephosY1068__FREE = kdephosY1068_f*species_50;  // dwdp[54]
    dflux_reaction_385_dkdephosY1068__FREE = kdephosY1068_f*species_51;  // dwdp[55]
    dflux_reaction_386_dkdephosY1068__FREE = kdephosY1068_f*species_55;  // dwdp[56]
    dflux_reaction_387_dkdephosY1068__FREE = kdephosY1068_f*species_57;  // dwdp[57]
    dflux_reaction_388_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_45;  // dwdp[58]
    dflux_reaction_389_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_47;  // dwdp[59]
    dflux_reaction_390_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_50;  // dwdp[60]
    dflux_reaction_391_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_52;  // dwdp[61]
    dflux_reaction_392_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_56;  // dwdp[62]
    dflux_reaction_393_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_58;  // dwdp[63]
    dflux_reaction_532_dkdephosY1068__FREE = kdephosY1068_f*species_61;  // dwdp[64]
    dflux_reaction_533_dkdephosY1068__FREE = kdephosY1068_f*species_62;  // dwdp[65]
    dflux_reaction_534_dkdephosY1068__FREE = kdephosY1068_f*species_66;  // dwdp[66]
    dflux_reaction_535_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_61;  // dwdp[67]
    dflux_reaction_536_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_63;  // dwdp[68]
    dflux_reaction_537_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_67;  // dwdp[69]
    dflux_reaction_602_dkdephosY1068__FREE = kdephosY1068_f*species_71;  // dwdp[70]
    dflux_reaction_603_dkdephosY1068__FREE = kdephosY1173_f*ratio_kdephosY1173*species_72;  // dwdp[71]
    dflux_reaction_8_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_6;  // dwdp[72]
    dflux_reaction_11_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_7;  // dwdp[73]
    dflux_reaction_22_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_8;  // dwdp[74]
    dflux_reaction_23_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_9;  // dwdp[75]
    dflux_reaction_28_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_11;  // dwdp[76]
    dflux_reaction_29_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_12;  // dwdp[77]
    dflux_reaction_34_dkdephosYN__FREE = species_10;  // dwdp[78]
    dflux_reaction_35_dkdephosYN__FREE = species_13;  // dwdp[79]
    dflux_reaction_72_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_17;  // dwdp[80]
    dflux_reaction_73_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_23;  // dwdp[81]
    dflux_reaction_74_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_25;  // dwdp[82]
    dflux_reaction_79_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_20;  // dwdp[83]
    dflux_reaction_80_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_24;  // dwdp[84]
    dflux_reaction_81_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_26;  // dwdp[85]
    dflux_reaction_92_dkdephosYN__FREE = species_16;  // dwdp[86]
    dflux_reaction_93_dkdephosYN__FREE = species_18;  // dwdp[87]
    dflux_reaction_94_dkdephosYN__FREE = species_19;  // dwdp[88]
    dflux_reaction_95_dkdephosYN__FREE = species_21;  // dwdp[89]
    dflux_reaction_96_dkdephosYN__FREE = species_22;  // dwdp[90]
    dflux_reaction_196_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_36;  // dwdp[91]
    dflux_reaction_197_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_37;  // dwdp[92]
    dflux_reaction_200_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_41;  // dwdp[93]
    dflux_reaction_201_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_42;  // dwdp[94]
    dflux_reaction_216_dkdephosYN__FREE = species_29;  // dwdp[95]
    dflux_reaction_217_dkdephosYN__FREE = species_31;  // dwdp[96]
    dflux_reaction_218_dkdephosYN__FREE = species_32;  // dwdp[97]
    dflux_reaction_219_dkdephosYN__FREE = species_35;  // dwdp[98]
    dflux_reaction_220_dkdephosYN__FREE = species_38;  // dwdp[99]
    dflux_reaction_221_dkdephosYN__FREE = species_39;  // dwdp[100]
    dflux_reaction_222_dkdephosYN__FREE = species_40;  // dwdp[101]
    dflux_reaction_223_dkdephosYN__FREE = species_43;  // dwdp[102]
    dflux_reaction_224_dkdephosYN__FREE = species_44;  // dwdp[103]
    dflux_reaction_380_dkdephosYN__FREE = kphos_f*ratio_kpkd_YN__FREE*species_59;  // dwdp[104]
    dflux_reaction_381_dkdephosYN__FREE = kphos_f*ratio_kphos_receiver*ratio_kpkd_YN__FREE*species_60;  // dwdp[105]
    dflux_reaction_394_dkdephosYN__FREE = species_46;  // dwdp[106]
    dflux_reaction_395_dkdephosYN__FREE = species_47;  // dwdp[107]
    dflux_reaction_396_dkdephosYN__FREE = species_50;  // dwdp[108]
    dflux_reaction_397_dkdephosYN__FREE = species_53;  // dwdp[109]
    dflux_reaction_398_dkdephosYN__FREE = species_54;  // dwdp[110]
    dflux_reaction_399_dkdephosYN__FREE = species_55;  // dwdp[111]
    dflux_reaction_400_dkdephosYN__FREE = species_56;  // dwdp[112]
    dflux_reaction_401_dkdephosYN__FREE = species_57;  // dwdp[113]
    dflux_reaction_402_dkdephosYN__FREE = species_58;  // dwdp[114]
    dflux_reaction_538_dkdephosYN__FREE = species_61;  // dwdp[115]
    dflux_reaction_539_dkdephosYN__FREE = species_64;  // dwdp[116]
    dflux_reaction_540_dkdephosYN__FREE = species_65;  // dwdp[117]
    dflux_reaction_541_dkdephosYN__FREE = species_66;  // dwdp[118]
    dflux_reaction_542_dkdephosYN__FREE = species_67;  // dwdp[119]
    dflux_reaction_543_dkdephosYN__FREE = species_69;  // dwdp[120]
    dflux_reaction_544_dkdephosYN__FREE = species_70;  // dwdp[121]
    dflux_reaction_604_dkdephosYN__FREE = species_71;  // dwdp[122]
    dflux_reaction_605_dkdephosYN__FREE = species_72;  // dwdp[123]
    dflux_reaction_606_dkdephosYN__FREE = species_74;  // dwdp[124]
    dflux_reaction_616_dkdephosYN__FREE = species_75;  // dwdp[125]
    dflux_reaction_6_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*species_6;  // dwdp[126]
    dflux_reaction_7_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*species_6;  // dwdp[127]
    dflux_reaction_9_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphos_receiver*species_7;  // dwdp[128]
    dflux_reaction_10_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*ratio_kphos_receiver*species_7;  // dwdp[129]
    dflux_reaction_18_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*species_9;  // dwdp[130]
    dflux_reaction_19_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*species_10;  // dwdp[131]
    dflux_reaction_20_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*species_8;  // dwdp[132]
    dflux_reaction_21_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*species_10;  // dwdp[133]
    dflux_reaction_24_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphos_receiver*species_12;  // dwdp[134]
    dflux_reaction_25_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphos_receiver*species_13;  // dwdp[135]
    dflux_reaction_26_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*ratio_kphos_receiver*species_11;  // dwdp[136]
    dflux_reaction_27_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*ratio_kphos_receiver*species_13;  // dwdp[137]
    dflux_reaction_68_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*species_19;  // dwdp[138]
    dflux_reaction_69_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*species_25;  // dwdp[139]
    dflux_reaction_70_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*species_18;  // dwdp[140]
    dflux_reaction_71_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*species_23;  // dwdp[141]
    dflux_reaction_75_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphos_receiver*species_22;  // dwdp[142]
    dflux_reaction_76_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphos_receiver*species_26;  // dwdp[143]
    dflux_reaction_77_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*ratio_kphos_receiver*species_21;  // dwdp[144]
    dflux_reaction_78_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*ratio_kphos_receiver*species_24;  // dwdp[145]
    dflux_reaction_194_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*species_39;  // dwdp[146]
    dflux_reaction_195_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*species_38;  // dwdp[147]
    dflux_reaction_198_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphos_receiver*species_44;  // dwdp[148]
    dflux_reaction_199_dratio_kpkd_Y1068__FREE = kdephosY1068__FREE*kphos_f*ratio_kphosY1173*ratio_kphos_receiver*species_43;  // dwdp[149]
    dflux_reaction_8_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_6;  // dwdp[150]
    dflux_reaction_11_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_7;  // dwdp[151]
    dflux_reaction_22_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_8;  // dwdp[152]
    dflux_reaction_23_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_9;  // dwdp[153]
    dflux_reaction_28_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_11;  // dwdp[154]
    dflux_reaction_29_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_12;  // dwdp[155]
    dflux_reaction_72_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_17;  // dwdp[156]
    dflux_reaction_73_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_23;  // dwdp[157]
    dflux_reaction_74_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_25;  // dwdp[158]
    dflux_reaction_79_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_20;  // dwdp[159]
    dflux_reaction_80_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_24;  // dwdp[160]
    dflux_reaction_81_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_26;  // dwdp[161]
    dflux_reaction_196_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_36;  // dwdp[162]
    dflux_reaction_197_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_37;  // dwdp[163]
    dflux_reaction_200_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_41;  // dwdp[164]
    dflux_reaction_201_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_42;  // dwdp[165]
    dflux_reaction_380_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*species_59;  // dwdp[166]
    dflux_reaction_381_dratio_kpkd_YN__FREE = kdephosYN__FREE*kphos_f*ratio_kphos_receiver*species_60;  // dwdp[167]
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
