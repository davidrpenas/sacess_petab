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
#include "xdot.h"
#include "xdot_old.h"
#include "sx.h"
#include "stau.h"

namespace amici {
namespace model_Smith_BMCSystBiol2013 {

void deltasx_Smith_BMCSystBiol2013(realtype *deltasx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const int ip, const int ie, const realtype *xdot, const realtype *xdot_old, const realtype *sx, const realtype *stau, const realtype *tcl){
    switch(ie) {
        case 0:
        case 3:
        case 4:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                case 19:
                case 20:
                case 21:
                case 22:
                case 23:
                case 24:
                case 25:
                case 26:
                case 27:
                    deltasx[1] = -sx1;
                    break;
            }
            break;
        case 1:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                case 19:
                case 20:
                case 21:
                case 22:
                case 23:
                case 24:
                case 25:
                case 26:
                case 27:
                    deltasx[0] = -sx11*(dNULLdt - xdot_old0)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[1] = -sx11*(dInsdt - xdot_old1)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[2] = -sx11*(dInRdt - xdot_old2)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[3] = -sx11*(dIns_InRdt - xdot_old3)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[4] = -sx11*(dIns_InR_Pdt - xdot_old4)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[5] = -sx11*(dIns_2_InR_Pdt - xdot_old5)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[6] = -sx11*(dcytoplasm_InRdt - xdot_old6)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[7] = -sx11*(dcytoplasm_Ins_2_InR_Pdt - xdot_old7)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[8] = -sx11*(dcytoplasm_Ins_InR_Pdt - xdot_old8)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[9] = -sx11*(dPTP1Bdt - xdot_old9)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[10] = -sx11*(dIRS1_TyrP_PI3Kdt - xdot_old10)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[11] = -sx11*(dPI345P3dt - xdot_old11)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[12] = -sx11*(dPIP2dt - xdot_old12)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[13] = -sx11*(dAktdt - xdot_old13)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[14] = -sx11*(dAkt_P2dt - xdot_old14)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[15] = -sx11*(dPKCdt - xdot_old15)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[16] = -sx11*(dPKC_Pdt - xdot_old16)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[17] = -sx11*(dAS160dt - xdot_old17)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[18] = -sx11*(dAS160_Pdt - xdot_old18)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[19] = -sx11*(dcytoplasm_GLUT4dt - xdot_old19)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[20] = -sx11*(dcellsurface_GLUT4dt - xdot_old20)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[21] = -sx11*(dPTENdt - xdot_old21)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[22] = -sx11*(dPP2Adt - xdot_old22)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[23] = -sx11*(dextracellular_ROSdt - xdot_old23)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[24] = -sx11*(dPTP1B_oxdt - xdot_old24)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[25] = -sx11*(dPTEN_oxdt - xdot_old25)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[26] = -sx11*(dROSdt - xdot_old26)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[27] = -sx11*(dGSHdt - xdot_old27)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[28] = -sx11*(dGSSGdt - xdot_old28)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[29] = -sx11*(dcytoplasm_SOD2dt - xdot_old29)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[30] = -sx11*(dNOX_inactdt - xdot_old30)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[31] = -sx11*(dNOXdt - xdot_old31)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[32] = -sx11*(dNOX_deactdt - xdot_old32)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[33] = -sx11*(dMtdt - xdot_old33)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[34] = -sx11*(dIRS1dt - xdot_old34)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[35] = -sx11*(dIRS1_TyrPdt - xdot_old35)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[36] = -sx11*(dIRS1_PolySerPdt - xdot_old36)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[37] = -sx11*(dIRS1_TyrP_PolySerPdt - xdot_old37)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[38] = -sx11*(dPI3Kdt - xdot_old38)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[39] = -sx11*(dJNK_Pdt - xdot_old39)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[40] = -sx11*(dIKK_Pdt - xdot_old40)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[41] = -sx11*(dJNKdt - xdot_old41)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[42] = -sx11*(dIKKdt - xdot_old42)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[43] = -sx11*(dDUSPdt - xdot_old43)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[44] = -sx11*(dDUSP_oxdt - xdot_old44)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[45] = -sx11*(dnulldt - xdot_old45)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[46] = -sx11*(ddegr_Foxo1dt - xdot_old46)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[47] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb0dt - xdot_old47)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[48] = -sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe0_pUb0dt - xdot_old48)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[49] = -sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe0_pUb0dt - xdot_old49)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[50] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb1dt - xdot_old50)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[51] = -sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe0_pUb1dt - xdot_old51)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[52] = -sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe0_pUb1dt - xdot_old52)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[53] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb0dt - xdot_old53)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[54] = -sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe1_pUb0dt - xdot_old54)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[55] = -sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe1_pUb0dt - xdot_old55)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[56] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb1dt - xdot_old56)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[57] = -sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe1_pUb1dt - xdot_old57)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[58] = -sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe1_pUb1dt - xdot_old58)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[59] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb0dt - xdot_old59)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[60] = -sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe0_pUb0dt - xdot_old60)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[61] = -sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe0_pUb0dt - xdot_old61)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[62] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb1dt - xdot_old62)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[63] = -sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe0_pUb1dt - xdot_old63)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[64] = -sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe0_pUb1dt - xdot_old64)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[65] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb0dt - xdot_old65)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[66] = -sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe1_pUb0dt - xdot_old66)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[67] = -sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe1_pUb0dt - xdot_old67)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[68] = -sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb1dt - xdot_old68)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[69] = -sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe1_pUb1dt - xdot_old69)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[70] = -sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe1_pUb1dt - xdot_old70)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[71] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb0dt - xdot_old71)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[72] = -sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe0_pUb0dt - xdot_old72)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[73] = -sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe0_pUb0dt - xdot_old73)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[74] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb1dt - xdot_old74)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[75] = -sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe0_pUb1dt - xdot_old75)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[76] = -sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe0_pUb1dt - xdot_old76)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[77] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb0dt - xdot_old77)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[78] = -sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe1_pUb0dt - xdot_old78)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[79] = -sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe1_pUb0dt - xdot_old79)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[80] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb1dt - xdot_old80)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[81] = -sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe1_pUb1dt - xdot_old81)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[82] = -sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe1_pUb1dt - xdot_old82)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[83] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb0dt - xdot_old83)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[84] = -sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe0_pUb0dt - xdot_old84)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[85] = -sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe0_pUb0dt - xdot_old85)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[86] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb1dt - xdot_old86)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[87] = -sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe0_pUb1dt - xdot_old87)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[88] = -sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe0_pUb1dt - xdot_old88)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[89] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb0dt - xdot_old89)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[90] = -sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe1_pUb0dt - xdot_old90)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[91] = -sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe1_pUb0dt - xdot_old91)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[92] = -sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb1dt - xdot_old92)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[93] = -sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe1_pUb1dt - xdot_old93)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[94] = -sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe1_pUb1dt - xdot_old94)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[95] = -sx11*(dnucleus_RNA_InRdt - xdot_old95)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[96] = -sx11*(dcytoplasm_RNA_InRdt - xdot_old96)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[97] = -sx11*(dnucleus_RNA_SOD2dt - xdot_old97)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[98] = -sx11*(dcytoplasm_RNA_SOD2dt - xdot_old98)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[99] = -sx11*(dE2F1dt - xdot_old99)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[100] = -sx11*(dSGKdt - xdot_old100)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[101] = -sx11*(dCDK2dt - xdot_old101)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[102] = -sx11*(dAMPKdt - xdot_old102)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[103] = -sx11*(dCBPP300dt - xdot_old103)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[104] = -sx11*(dSIRT1dt - xdot_old104)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[105] = -sx11*(dE3dt - xdot_old105)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[106] = -sx11*(dUSP7dt - xdot_old106)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[107] = -sx11*(dSCFdt - xdot_old107)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[108] = -sx11*(dProteasomedt - xdot_old108)/(PI345P3*(PTEN*kminus9 + kminus9_basal) - PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    break;
            }
            break;
        case 2:
            switch(ip) {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                case 19:
                case 20:
                case 21:
                case 22:
                case 23:
                case 24:
                case 25:
                case 26:
                case 27:
                    deltasx[0] = sx11*(dNULLdt - xdot_old0)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[1] = sx11*(dInsdt - xdot_old1)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[2] = sx11*(dInRdt - xdot_old2)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[3] = sx11*(dIns_InRdt - xdot_old3)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[4] = sx11*(dIns_InR_Pdt - xdot_old4)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[5] = sx11*(dIns_2_InR_Pdt - xdot_old5)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[6] = sx11*(dcytoplasm_InRdt - xdot_old6)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[7] = sx11*(dcytoplasm_Ins_2_InR_Pdt - xdot_old7)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[8] = sx11*(dcytoplasm_Ins_InR_Pdt - xdot_old8)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[9] = sx11*(dPTP1Bdt - xdot_old9)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[10] = sx11*(dIRS1_TyrP_PI3Kdt - xdot_old10)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[11] = sx11*(dPI345P3dt - xdot_old11)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[12] = sx11*(dPIP2dt - xdot_old12)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[13] = sx11*(dAktdt - xdot_old13)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[14] = sx11*(dAkt_P2dt - xdot_old14)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[15] = sx11*(dPKCdt - xdot_old15)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[16] = sx11*(dPKC_Pdt - xdot_old16)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[17] = sx11*(dAS160dt - xdot_old17)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[18] = sx11*(dAS160_Pdt - xdot_old18)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[19] = sx11*(dcytoplasm_GLUT4dt - xdot_old19)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[20] = sx11*(dcellsurface_GLUT4dt - xdot_old20)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[21] = sx11*(dPTENdt - xdot_old21)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[22] = sx11*(dPP2Adt - xdot_old22)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[23] = sx11*(dextracellular_ROSdt - xdot_old23)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[24] = sx11*(dPTP1B_oxdt - xdot_old24)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[25] = sx11*(dPTEN_oxdt - xdot_old25)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[26] = sx11*(dROSdt - xdot_old26)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[27] = sx11*(dGSHdt - xdot_old27)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[28] = sx11*(dGSSGdt - xdot_old28)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[29] = sx11*(dcytoplasm_SOD2dt - xdot_old29)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[30] = sx11*(dNOX_inactdt - xdot_old30)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[31] = sx11*(dNOXdt - xdot_old31)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[32] = sx11*(dNOX_deactdt - xdot_old32)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[33] = sx11*(dMtdt - xdot_old33)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[34] = sx11*(dIRS1dt - xdot_old34)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[35] = sx11*(dIRS1_TyrPdt - xdot_old35)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[36] = sx11*(dIRS1_PolySerPdt - xdot_old36)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[37] = sx11*(dIRS1_TyrP_PolySerPdt - xdot_old37)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[38] = sx11*(dPI3Kdt - xdot_old38)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[39] = sx11*(dJNK_Pdt - xdot_old39)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[40] = sx11*(dIKK_Pdt - xdot_old40)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[41] = sx11*(dJNKdt - xdot_old41)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[42] = sx11*(dIKKdt - xdot_old42)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[43] = sx11*(dDUSPdt - xdot_old43)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[44] = sx11*(dDUSP_oxdt - xdot_old44)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[45] = sx11*(dnulldt - xdot_old45)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[46] = sx11*(ddegr_Foxo1dt - xdot_old46)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[47] = sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb0dt - xdot_old47)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[48] = sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe0_pUb0dt - xdot_old48)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[49] = sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe0_pUb0dt - xdot_old49)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[50] = sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe0_pUb1dt - xdot_old50)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[51] = sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe0_pUb1dt - xdot_old51)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[52] = sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe0_pUb1dt - xdot_old52)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[53] = sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb0dt - xdot_old53)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[54] = sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe1_pUb0dt - xdot_old54)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[55] = sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe1_pUb0dt - xdot_old55)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[56] = sx11*(dcytoplasm_Foxo1_Pa0_Pd0_Pe1_pUb1dt - xdot_old56)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[57] = sx11*(dnucleus_Foxo1_Pa0_Pd0_Pe1_pUb1dt - xdot_old57)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[58] = sx11*(ddnabound_Foxo1_Pa0_Pd0_Pe1_pUb1dt - xdot_old58)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[59] = sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb0dt - xdot_old59)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[60] = sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe0_pUb0dt - xdot_old60)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[61] = sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe0_pUb0dt - xdot_old61)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[62] = sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe0_pUb1dt - xdot_old62)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[63] = sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe0_pUb1dt - xdot_old63)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[64] = sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe0_pUb1dt - xdot_old64)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[65] = sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb0dt - xdot_old65)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[66] = sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe1_pUb0dt - xdot_old66)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[67] = sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe1_pUb0dt - xdot_old67)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[68] = sx11*(dcytoplasm_Foxo1_Pa0_Pd1_Pe1_pUb1dt - xdot_old68)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[69] = sx11*(dnucleus_Foxo1_Pa0_Pd1_Pe1_pUb1dt - xdot_old69)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[70] = sx11*(ddnabound_Foxo1_Pa0_Pd1_Pe1_pUb1dt - xdot_old70)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[71] = sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb0dt - xdot_old71)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[72] = sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe0_pUb0dt - xdot_old72)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[73] = sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe0_pUb0dt - xdot_old73)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[74] = sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe0_pUb1dt - xdot_old74)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[75] = sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe0_pUb1dt - xdot_old75)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[76] = sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe0_pUb1dt - xdot_old76)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[77] = sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb0dt - xdot_old77)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[78] = sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe1_pUb0dt - xdot_old78)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[79] = sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe1_pUb0dt - xdot_old79)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[80] = sx11*(dcytoplasm_Foxo1_Pa1_Pd0_Pe1_pUb1dt - xdot_old80)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[81] = sx11*(dnucleus_Foxo1_Pa1_Pd0_Pe1_pUb1dt - xdot_old81)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[82] = sx11*(ddnabound_Foxo1_Pa1_Pd0_Pe1_pUb1dt - xdot_old82)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[83] = sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb0dt - xdot_old83)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[84] = sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe0_pUb0dt - xdot_old84)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[85] = sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe0_pUb0dt - xdot_old85)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[86] = sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe0_pUb1dt - xdot_old86)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[87] = sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe0_pUb1dt - xdot_old87)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[88] = sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe0_pUb1dt - xdot_old88)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[89] = sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb0dt - xdot_old89)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[90] = sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe1_pUb0dt - xdot_old90)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[91] = sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe1_pUb0dt - xdot_old91)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[92] = sx11*(dcytoplasm_Foxo1_Pa1_Pd1_Pe1_pUb1dt - xdot_old92)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[93] = sx11*(dnucleus_Foxo1_Pa1_Pd1_Pe1_pUb1dt - xdot_old93)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[94] = sx11*(ddnabound_Foxo1_Pa1_Pd1_Pe1_pUb1dt - xdot_old94)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[95] = sx11*(dnucleus_RNA_InRdt - xdot_old95)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[96] = sx11*(dcytoplasm_RNA_InRdt - xdot_old96)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[97] = sx11*(dnucleus_RNA_SOD2dt - xdot_old97)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[98] = sx11*(dcytoplasm_RNA_SOD2dt - xdot_old98)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[99] = sx11*(dE2F1dt - xdot_old99)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[100] = sx11*(dSGKdt - xdot_old100)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[101] = sx11*(dCDK2dt - xdot_old101)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[102] = sx11*(dAMPKdt - xdot_old102)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[103] = sx11*(dCBPP300dt - xdot_old103)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[104] = sx11*(dSIRT1dt - xdot_old104)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[105] = sx11*(dE3dt - xdot_old105)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[106] = sx11*(dUSP7dt - xdot_old106)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[107] = sx11*(dSCFdt - xdot_old107)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    deltasx[108] = sx11*(dProteasomedt - xdot_old108)/(-PI345P3*(PTEN*kminus9 + kminus9_basal) + PIP2*(IRS1_TyrP_PI3K*k9 + k9_basal));
                    break;
            }
            break;
    }
}

} // namespace model_Smith_BMCSystBiol2013
} // namespace amici
