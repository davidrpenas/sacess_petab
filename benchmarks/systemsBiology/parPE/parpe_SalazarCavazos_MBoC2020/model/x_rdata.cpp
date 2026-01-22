#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void x_rdata_SalazarCavazos_MBoC2020(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = species_1;
    x_rdata[1] = species_2;
    x_rdata[2] = species_3;
    x_rdata[3] = species_4;
    x_rdata[4] = species_5;
    x_rdata[5] = species_6;
    x_rdata[6] = species_7;
    x_rdata[7] = species_8;
    x_rdata[8] = species_9;
    x_rdata[9] = species_10;
    x_rdata[10] = species_11;
    x_rdata[11] = species_12;
    x_rdata[12] = species_13;
    x_rdata[13] = species_14;
    x_rdata[14] = species_15;
    x_rdata[15] = species_16;
    x_rdata[16] = species_17;
    x_rdata[17] = species_18;
    x_rdata[18] = species_19;
    x_rdata[19] = species_20;
    x_rdata[20] = species_21;
    x_rdata[21] = species_22;
    x_rdata[22] = species_23;
    x_rdata[23] = species_24;
    x_rdata[24] = species_25;
    x_rdata[25] = species_26;
    x_rdata[26] = species_27;
    x_rdata[27] = species_28;
    x_rdata[28] = species_29;
    x_rdata[29] = species_30;
    x_rdata[30] = species_31;
    x_rdata[31] = species_32;
    x_rdata[32] = species_33;
    x_rdata[33] = species_34;
    x_rdata[34] = species_35;
    x_rdata[35] = species_36;
    x_rdata[36] = species_37;
    x_rdata[37] = species_38;
    x_rdata[38] = species_39;
    x_rdata[39] = species_40;
    x_rdata[40] = species_41;
    x_rdata[41] = species_42;
    x_rdata[42] = species_43;
    x_rdata[43] = species_44;
    x_rdata[44] = species_45;
    x_rdata[45] = species_46;
    x_rdata[46] = species_47;
    x_rdata[47] = species_48;
    x_rdata[48] = species_49;
    x_rdata[49] = species_50;
    x_rdata[50] = species_51;
    x_rdata[51] = species_52;
    x_rdata[52] = species_53;
    x_rdata[53] = species_54;
    x_rdata[54] = species_55;
    x_rdata[55] = species_56;
    x_rdata[56] = species_57;
    x_rdata[57] = species_58;
    x_rdata[58] = species_59;
    x_rdata[59] = species_60;
    x_rdata[60] = species_61;
    x_rdata[61] = species_62;
    x_rdata[62] = species_63;
    x_rdata[63] = species_64;
    x_rdata[64] = species_65;
    x_rdata[65] = species_66;
    x_rdata[66] = species_67;
    x_rdata[67] = species_68;
    x_rdata[68] = species_69;
    x_rdata[69] = species_70;
    x_rdata[70] = species_71;
    x_rdata[71] = species_72;
    x_rdata[72] = species_73;
    x_rdata[73] = species_74;
    x_rdata[74] = species_75;
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
