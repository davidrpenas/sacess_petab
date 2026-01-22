#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x_rdata.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void x_solver_Smith_BMCSystBiol2013(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = NULL;
    x_solver[1] = Ins;
    x_solver[2] = InR;
    x_solver[3] = Ins_InR;
    x_solver[4] = Ins_InR_P;
    x_solver[5] = Ins_2_InR_P;
    x_solver[6] = cytoplasm_InR;
    x_solver[7] = cytoplasm_Ins_2_InR_P;
    x_solver[8] = cytoplasm_Ins_InR_P;
    x_solver[9] = PTP1B;
    x_solver[10] = IRS1_TyrP_PI3K;
    x_solver[11] = PI345P3;
    x_solver[12] = PIP2;
    x_solver[13] = Akt;
    x_solver[14] = Akt_P2;
    x_solver[15] = PKC;
    x_solver[16] = PKC_P;
    x_solver[17] = AS160;
    x_solver[18] = AS160_P;
    x_solver[19] = cytoplasm_GLUT4;
    x_solver[20] = cellsurface_GLUT4;
    x_solver[21] = PTEN;
    x_solver[22] = PP2A;
    x_solver[23] = extracellular_ROS;
    x_solver[24] = PTP1B_ox;
    x_solver[25] = PTEN_ox;
    x_solver[26] = ROS;
    x_solver[27] = GSH;
    x_solver[28] = GSSG;
    x_solver[29] = cytoplasm_SOD2;
    x_solver[30] = NOX_inact;
    x_solver[31] = NOX;
    x_solver[32] = NOX_deact;
    x_solver[33] = Mt;
    x_solver[34] = IRS1;
    x_solver[35] = IRS1_TyrP;
    x_solver[36] = IRS1_PolySerP;
    x_solver[37] = IRS1_TyrP_PolySerP;
    x_solver[38] = PI3K;
    x_solver[39] = JNK_P;
    x_solver[40] = IKK_P;
    x_solver[41] = JNK;
    x_solver[42] = IKK;
    x_solver[43] = DUSP;
    x_solver[44] = DUSP_ox;
    x_solver[45] = null;
    x_solver[46] = degr_Foxo1;
    x_solver[47] = cytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb0;
    x_solver[48] = nucleus_Foxo1_Pa0_Pd0_Pe0_pUb0;
    x_solver[49] = dnabound_Foxo1_Pa0_Pd0_Pe0_pUb0;
    x_solver[50] = cytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb1;
    x_solver[51] = nucleus_Foxo1_Pa0_Pd0_Pe0_pUb1;
    x_solver[52] = dnabound_Foxo1_Pa0_Pd0_Pe0_pUb1;
    x_solver[53] = cytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb0;
    x_solver[54] = nucleus_Foxo1_Pa0_Pd0_Pe1_pUb0;
    x_solver[55] = dnabound_Foxo1_Pa0_Pd0_Pe1_pUb0;
    x_solver[56] = cytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb1;
    x_solver[57] = nucleus_Foxo1_Pa0_Pd0_Pe1_pUb1;
    x_solver[58] = dnabound_Foxo1_Pa0_Pd0_Pe1_pUb1;
    x_solver[59] = cytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb0;
    x_solver[60] = nucleus_Foxo1_Pa0_Pd1_Pe0_pUb0;
    x_solver[61] = dnabound_Foxo1_Pa0_Pd1_Pe0_pUb0;
    x_solver[62] = cytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb1;
    x_solver[63] = nucleus_Foxo1_Pa0_Pd1_Pe0_pUb1;
    x_solver[64] = dnabound_Foxo1_Pa0_Pd1_Pe0_pUb1;
    x_solver[65] = cytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb0;
    x_solver[66] = nucleus_Foxo1_Pa0_Pd1_Pe1_pUb0;
    x_solver[67] = dnabound_Foxo1_Pa0_Pd1_Pe1_pUb0;
    x_solver[68] = cytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb1;
    x_solver[69] = nucleus_Foxo1_Pa0_Pd1_Pe1_pUb1;
    x_solver[70] = dnabound_Foxo1_Pa0_Pd1_Pe1_pUb1;
    x_solver[71] = cytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb0;
    x_solver[72] = nucleus_Foxo1_Pa1_Pd0_Pe0_pUb0;
    x_solver[73] = dnabound_Foxo1_Pa1_Pd0_Pe0_pUb0;
    x_solver[74] = cytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb1;
    x_solver[75] = nucleus_Foxo1_Pa1_Pd0_Pe0_pUb1;
    x_solver[76] = dnabound_Foxo1_Pa1_Pd0_Pe0_pUb1;
    x_solver[77] = cytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb0;
    x_solver[78] = nucleus_Foxo1_Pa1_Pd0_Pe1_pUb0;
    x_solver[79] = dnabound_Foxo1_Pa1_Pd0_Pe1_pUb0;
    x_solver[80] = cytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb1;
    x_solver[81] = nucleus_Foxo1_Pa1_Pd0_Pe1_pUb1;
    x_solver[82] = dnabound_Foxo1_Pa1_Pd0_Pe1_pUb1;
    x_solver[83] = cytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb0;
    x_solver[84] = nucleus_Foxo1_Pa1_Pd1_Pe0_pUb0;
    x_solver[85] = dnabound_Foxo1_Pa1_Pd1_Pe0_pUb0;
    x_solver[86] = cytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb1;
    x_solver[87] = nucleus_Foxo1_Pa1_Pd1_Pe0_pUb1;
    x_solver[88] = dnabound_Foxo1_Pa1_Pd1_Pe0_pUb1;
    x_solver[89] = cytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb0;
    x_solver[90] = nucleus_Foxo1_Pa1_Pd1_Pe1_pUb0;
    x_solver[91] = dnabound_Foxo1_Pa1_Pd1_Pe1_pUb0;
    x_solver[92] = cytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb1;
    x_solver[93] = nucleus_Foxo1_Pa1_Pd1_Pe1_pUb1;
    x_solver[94] = dnabound_Foxo1_Pa1_Pd1_Pe1_pUb1;
    x_solver[95] = nucleus_RNA_InR;
    x_solver[96] = cytoplasm_RNA_InR;
    x_solver[97] = nucleus_RNA_SOD2;
    x_solver[98] = cytoplasm_RNA_SOD2;
    x_solver[99] = E2F1;
    x_solver[100] = SGK;
    x_solver[101] = CDK2;
    x_solver[102] = AMPK;
    x_solver[103] = CBPP300;
    x_solver[104] = SIRT1;
    x_solver[105] = E3;
    x_solver[106] = USP7;
    x_solver[107] = SCF;
    x_solver[108] = Proteasome;
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
