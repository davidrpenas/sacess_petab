#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_SalazarCavazos_MBoC2020 {

void x_solver_SalazarCavazos_MBoC2020(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = species_1;
    x_solver[1] = species_2;
    x_solver[2] = species_3;
    x_solver[3] = species_4;
    x_solver[4] = species_5;
    x_solver[5] = species_6;
    x_solver[6] = species_7;
    x_solver[7] = species_8;
    x_solver[8] = species_9;
    x_solver[9] = species_10;
    x_solver[10] = species_11;
    x_solver[11] = species_12;
    x_solver[12] = species_13;
    x_solver[13] = species_14;
    x_solver[14] = species_15;
    x_solver[15] = species_16;
    x_solver[16] = species_17;
    x_solver[17] = species_18;
    x_solver[18] = species_19;
    x_solver[19] = species_20;
    x_solver[20] = species_21;
    x_solver[21] = species_22;
    x_solver[22] = species_23;
    x_solver[23] = species_24;
    x_solver[24] = species_25;
    x_solver[25] = species_26;
    x_solver[26] = species_27;
    x_solver[27] = species_28;
    x_solver[28] = species_29;
    x_solver[29] = species_30;
    x_solver[30] = species_31;
    x_solver[31] = species_32;
    x_solver[32] = species_33;
    x_solver[33] = species_34;
    x_solver[34] = species_35;
    x_solver[35] = species_36;
    x_solver[36] = species_37;
    x_solver[37] = species_38;
    x_solver[38] = species_39;
    x_solver[39] = species_40;
    x_solver[40] = species_41;
    x_solver[41] = species_42;
    x_solver[42] = species_43;
    x_solver[43] = species_44;
    x_solver[44] = species_45;
    x_solver[45] = species_46;
    x_solver[46] = species_47;
    x_solver[47] = species_48;
    x_solver[48] = species_49;
    x_solver[49] = species_50;
    x_solver[50] = species_51;
    x_solver[51] = species_52;
    x_solver[52] = species_53;
    x_solver[53] = species_54;
    x_solver[54] = species_55;
    x_solver[55] = species_56;
    x_solver[56] = species_57;
    x_solver[57] = species_58;
    x_solver[58] = species_59;
    x_solver[59] = species_60;
    x_solver[60] = species_61;
    x_solver[61] = species_62;
    x_solver[62] = species_63;
    x_solver[63] = species_64;
    x_solver[64] = species_65;
    x_solver[65] = species_66;
    x_solver[66] = species_67;
    x_solver[67] = species_68;
    x_solver[68] = species_69;
    x_solver[69] = species_70;
    x_solver[70] = species_71;
    x_solver[71] = species_72;
    x_solver[72] = species_73;
    x_solver[73] = species_74;
    x_solver[74] = species_75;
}

} // namespace model_SalazarCavazos_MBoC2020
} // namespace amici
