#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "x.h"
#include "p.h"
#include "k.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void x_rdata_Smith_BMCSystBiol2013(realtype *x_rdata, const realtype *x, const realtype *tcl, const realtype *p, const realtype *k){
    x_rdata[0] = NULL;
    x_rdata[1] = Ins;
    x_rdata[2] = InR;
    x_rdata[3] = Ins_InR;
    x_rdata[4] = Ins_InR_P;
    x_rdata[5] = Ins_2_InR_P;
    x_rdata[6] = cytoplasm_InR;
    x_rdata[7] = cytoplasm_Ins_2_InR_P;
    x_rdata[8] = cytoplasm_Ins_InR_P;
    x_rdata[9] = PTP1B;
    x_rdata[10] = IRS1_TyrP_PI3K;
    x_rdata[11] = PI345P3;
    x_rdata[12] = PIP2;
    x_rdata[13] = Akt;
    x_rdata[14] = Akt_P2;
    x_rdata[15] = PKC;
    x_rdata[16] = PKC_P;
    x_rdata[17] = AS160;
    x_rdata[18] = AS160_P;
    x_rdata[19] = cytoplasm_GLUT4;
    x_rdata[20] = cellsurface_GLUT4;
    x_rdata[21] = PTEN;
    x_rdata[22] = PP2A;
    x_rdata[23] = extracellular_ROS;
    x_rdata[24] = PTP1B_ox;
    x_rdata[25] = PTEN_ox;
    x_rdata[26] = ROS;
    x_rdata[27] = GSH;
    x_rdata[28] = GSSG;
    x_rdata[29] = cytoplasm_SOD2;
    x_rdata[30] = NOX_inact;
    x_rdata[31] = NOX;
    x_rdata[32] = NOX_deact;
    x_rdata[33] = Mt;
    x_rdata[34] = IRS1;
    x_rdata[35] = IRS1_TyrP;
    x_rdata[36] = IRS1_PolySerP;
    x_rdata[37] = IRS1_TyrP_PolySerP;
    x_rdata[38] = PI3K;
    x_rdata[39] = JNK_P;
    x_rdata[40] = IKK_P;
    x_rdata[41] = JNK;
    x_rdata[42] = IKK;
    x_rdata[43] = DUSP;
    x_rdata[44] = DUSP_ox;
    x_rdata[45] = null;
    x_rdata[46] = degr_Foxo1;
    x_rdata[47] = cytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb0;
    x_rdata[48] = nucleus_Foxo1_Pa0_Pd0_Pe0_pUb0;
    x_rdata[49] = dnabound_Foxo1_Pa0_Pd0_Pe0_pUb0;
    x_rdata[50] = cytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb1;
    x_rdata[51] = nucleus_Foxo1_Pa0_Pd0_Pe0_pUb1;
    x_rdata[52] = dnabound_Foxo1_Pa0_Pd0_Pe0_pUb1;
    x_rdata[53] = cytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb0;
    x_rdata[54] = nucleus_Foxo1_Pa0_Pd0_Pe1_pUb0;
    x_rdata[55] = dnabound_Foxo1_Pa0_Pd0_Pe1_pUb0;
    x_rdata[56] = cytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb1;
    x_rdata[57] = nucleus_Foxo1_Pa0_Pd0_Pe1_pUb1;
    x_rdata[58] = dnabound_Foxo1_Pa0_Pd0_Pe1_pUb1;
    x_rdata[59] = cytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb0;
    x_rdata[60] = nucleus_Foxo1_Pa0_Pd1_Pe0_pUb0;
    x_rdata[61] = dnabound_Foxo1_Pa0_Pd1_Pe0_pUb0;
    x_rdata[62] = cytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb1;
    x_rdata[63] = nucleus_Foxo1_Pa0_Pd1_Pe0_pUb1;
    x_rdata[64] = dnabound_Foxo1_Pa0_Pd1_Pe0_pUb1;
    x_rdata[65] = cytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb0;
    x_rdata[66] = nucleus_Foxo1_Pa0_Pd1_Pe1_pUb0;
    x_rdata[67] = dnabound_Foxo1_Pa0_Pd1_Pe1_pUb0;
    x_rdata[68] = cytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb1;
    x_rdata[69] = nucleus_Foxo1_Pa0_Pd1_Pe1_pUb1;
    x_rdata[70] = dnabound_Foxo1_Pa0_Pd1_Pe1_pUb1;
    x_rdata[71] = cytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb0;
    x_rdata[72] = nucleus_Foxo1_Pa1_Pd0_Pe0_pUb0;
    x_rdata[73] = dnabound_Foxo1_Pa1_Pd0_Pe0_pUb0;
    x_rdata[74] = cytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb1;
    x_rdata[75] = nucleus_Foxo1_Pa1_Pd0_Pe0_pUb1;
    x_rdata[76] = dnabound_Foxo1_Pa1_Pd0_Pe0_pUb1;
    x_rdata[77] = cytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb0;
    x_rdata[78] = nucleus_Foxo1_Pa1_Pd0_Pe1_pUb0;
    x_rdata[79] = dnabound_Foxo1_Pa1_Pd0_Pe1_pUb0;
    x_rdata[80] = cytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb1;
    x_rdata[81] = nucleus_Foxo1_Pa1_Pd0_Pe1_pUb1;
    x_rdata[82] = dnabound_Foxo1_Pa1_Pd0_Pe1_pUb1;
    x_rdata[83] = cytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb0;
    x_rdata[84] = nucleus_Foxo1_Pa1_Pd1_Pe0_pUb0;
    x_rdata[85] = dnabound_Foxo1_Pa1_Pd1_Pe0_pUb0;
    x_rdata[86] = cytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb1;
    x_rdata[87] = nucleus_Foxo1_Pa1_Pd1_Pe0_pUb1;
    x_rdata[88] = dnabound_Foxo1_Pa1_Pd1_Pe0_pUb1;
    x_rdata[89] = cytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb0;
    x_rdata[90] = nucleus_Foxo1_Pa1_Pd1_Pe1_pUb0;
    x_rdata[91] = dnabound_Foxo1_Pa1_Pd1_Pe1_pUb0;
    x_rdata[92] = cytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb1;
    x_rdata[93] = nucleus_Foxo1_Pa1_Pd1_Pe1_pUb1;
    x_rdata[94] = dnabound_Foxo1_Pa1_Pd1_Pe1_pUb1;
    x_rdata[95] = nucleus_RNA_InR;
    x_rdata[96] = cytoplasm_RNA_InR;
    x_rdata[97] = nucleus_RNA_SOD2;
    x_rdata[98] = cytoplasm_RNA_SOD2;
    x_rdata[99] = E2F1;
    x_rdata[100] = SGK;
    x_rdata[101] = CDK2;
    x_rdata[102] = AMPK;
    x_rdata[103] = CBPP300;
    x_rdata[104] = SIRT1;
    x_rdata[105] = E3;
    x_rdata[106] = USP7;
    x_rdata[107] = SCF;
    x_rdata[108] = Proteasome;
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
