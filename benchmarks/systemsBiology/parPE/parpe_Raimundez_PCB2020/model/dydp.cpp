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

namespace amici {
namespace model_Raimundez_PCB2020 {

void dydp_Raimundez_PCB2020(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 0:
            dydp[79] = 1/(KD_EGFR_CET*M_LN10);
            break;
        case 1:
            dydp[80] = 1/(KD_EGFR_EGF*M_LN10);
            break;
        case 3:
            dydp[81] = 1;
            break;
        case 4:
            dydp[82] = 1;
            break;
        case 5:
            dydp[83] = 1;
            break;
        case 6:
            dydp[84] = 1;
            break;
        case 7:
            dydp[85] = 1;
            break;
        case 8:
            dydp[86] = 1;
            break;
        case 9:
            dydp[87] = 1;
            break;
        case 10:
            dydp[88] = 1;
            break;
        case 11:
            dydp[89] = 1;
            break;
        case 12:
            dydp[90] = 1;
            break;
        case 13:
            dydp[91] = 1;
            break;
        case 14:
            dydp[92] = 1;
            break;
        case 15:
            dydp[93] = 1;
            break;
        case 16:
            dydp[94] = 1;
            break;
        case 17:
            dydp[95] = 1;
            break;
        case 33:
            dydp[96] = 1/(ki_AKT__MKN1*M_LN10);
            break;
        case 37:
            dydp[97] = 1/(kimp_pEGFR_EGF_2__MKN1*M_LN10);
            break;
        case 42:
            dydp[98] = 1/(xi_ka_PI3K_pMMET_2*M_LN10);
            break;
        case 43:
            dydp[99] = 1/(xi_ka_PI3K_pMMET_pEGFR*M_LN10);
            break;
        case 44:
            dydp[100] = 1/(xi_ka_RAS_pMMET_2*M_LN10);
            break;
        case 45:
            dydp[101] = 1/(xi_ka_RAS_pMMET_pEGFR*M_LN10);
            break;
        case 46:
            dydp[102] = 1/(xi_kdeg_pMMET_2_i*M_LN10);
            break;
        case 47:
            dydp[103] = 1/(xi_kdeg_pMMET_pEGFR_i*M_LN10);
            break;
        case 48:
            dydp[104] = 1/(xi_kdim_MMET*M_LN10);
            break;
        case 49:
            dydp[105] = 1/(xi_kdim_MMET_EGFR*M_LN10);
            break;
        case 50:
            dydp[106] = 1/(xi_kexp_pMMET_2_i*M_LN10);
            break;
        case 51:
            dydp[107] = 1/(xi_kexp_pMMET_pEGFR_i*M_LN10);
            break;
        case 52:
            dydp[108] = 1/(xi_ki_MPI3K*M_LN10);
            break;
        case 53:
            dydp[109] = 1/(xi_kimp_pMMET_2*M_LN10);
            break;
        case 54:
            dydp[110] = 1/(xi_kimp_pMMET_pEGFR*M_LN10);
            break;
        case 57:
            dydp[0] = EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 58:
            dydp[1] = EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 59:
            dydp[2] = EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 60:
            dydp[3] = EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 61:
            dydp[4] = EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 62:
            dydp[5] = EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 63:
            dydp[6] = EGFR + EGFR_CET + EGFR_EGF + 2*EGFR_EGF_2 + EGFR_MMET_METinh + MMET_EGFR + 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 64:
            dydp[7] = EGFR_MMET_METinh + MMET + 2*MMET_2 + MMET_EGFR + MMET_METinh + 2*MMET_METinh_2 + 2*MMET_MMET_METinh + 2*pMMET_2 + 2*pMMET_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 65:
            dydp[8] = EGFR_MMET_METinh + MMET + 2*MMET_2 + MMET_EGFR + MMET_METinh + 2*MMET_METinh_2 + 2*MMET_MMET_METinh + 2*pMMET_2 + 2*pMMET_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 66:
            dydp[9] = pAKT;
            break;
        case 67:
            dydp[10] = pAKT;
            break;
        case 68:
            dydp[11] = pAKT;
            break;
        case 69:
            dydp[12] = pAKT;
            break;
        case 70:
            dydp[13] = pAKT;
            break;
        case 71:
            dydp[14] = pAKT;
            break;
        case 72:
            dydp[15] = pAKT;
            break;
        case 73:
            dydp[16] = pAKT;
            break;
        case 74:
            dydp[17] = pAKT;
            break;
        case 75:
            dydp[18] = pAKT;
            break;
        case 76:
            dydp[19] = pAKT;
            break;
        case 77:
            dydp[20] = pAKT;
            break;
        case 78:
            dydp[21] = pAKT;
            break;
        case 79:
            dydp[22] = pAKT;
            break;
        case 80:
            dydp[23] = pAKT;
            break;
        case 81:
            dydp[24] = pAKT;
            break;
        case 82:
            dydp[25] = pAKT;
            break;
        case 83:
            dydp[26] = pAKT;
            break;
        case 84:
            dydp[27] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 85:
            dydp[28] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 86:
            dydp[29] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 87:
            dydp[30] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 88:
            dydp[31] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 89:
            dydp[32] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 90:
            dydp[33] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 91:
            dydp[34] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 92:
            dydp[35] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 93:
            dydp[36] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 94:
            dydp[37] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 95:
            dydp[38] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 96:
            dydp[39] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 97:
            dydp[40] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 98:
            dydp[41] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 99:
            dydp[42] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 100:
            dydp[43] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 101:
            dydp[44] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 102:
            dydp[45] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 103:
            dydp[46] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 104:
            dydp[47] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 105:
            dydp[48] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 106:
            dydp[49] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 107:
            dydp[50] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 108:
            dydp[51] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 109:
            dydp[52] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 110:
            dydp[53] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 111:
            dydp[54] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 112:
            dydp[55] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 113:
            dydp[56] = 2*pEGFR_EGF_2 + 2*pEGFR_EGF_2_i + pMMET_pEGFR + pMMET_pEGFR_i;
            break;
        case 114:
            dydp[57] = pMAPK;
            break;
        case 115:
            dydp[58] = pMAPK;
            break;
        case 116:
            dydp[59] = pMAPK;
            break;
        case 117:
            dydp[60] = pMAPK;
            break;
        case 118:
            dydp[61] = pMAPK;
            break;
        case 119:
            dydp[62] = pMAPK;
            break;
        case 120:
            dydp[63] = pMAPK;
            break;
        case 121:
            dydp[64] = pMAPK;
            break;
        case 122:
            dydp[65] = pMAPK;
            break;
        case 123:
            dydp[66] = pMAPK;
            break;
        case 124:
            dydp[67] = pMAPK;
            break;
        case 125:
            dydp[68] = pMAPK;
            break;
        case 126:
            dydp[69] = pMAPK;
            break;
        case 127:
            dydp[70] = pMAPK;
            break;
        case 128:
            dydp[71] = pMAPK;
            break;
        case 129:
            dydp[72] = pMAPK;
            break;
        case 130:
            dydp[73] = pMAPK;
            break;
        case 131:
            dydp[74] = pMAPK;
            break;
        case 132:
            dydp[75] = pMAPK;
            break;
        case 133:
            dydp[76] = pMAPK;
            break;
        case 134:
            dydp[77] = pMAPK;
            break;
        case 135:
            dydp[78] = pMAPK;
            break;
    }
}

} // namespace model_Raimundez_PCB2020
} // namespace amici
