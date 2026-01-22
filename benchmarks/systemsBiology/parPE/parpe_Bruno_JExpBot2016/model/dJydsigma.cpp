#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"
#include "my.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void dJydsigma_Bruno_JExpBot2016(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_ob10__0_043179__model1_data4 - 1.0*std::pow(-mob10__0_043179__model1_data4 + ob10__0_043179__model1_data4, 2)/std::pow(sigma_ob10__0_043179__model1_data4, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_ob10__0_046576__model1_data4 - 1.0*std::pow(-mob10__0_046576__model1_data4 + ob10__0_046576__model1_data4, 2)/std::pow(sigma_ob10__0_046576__model1_data4, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_ob10__0_048771__model1_data4 - 1.0*std::pow(-mob10__0_048771__model1_data4 + ob10__0_048771__model1_data4, 2)/std::pow(sigma_ob10__0_048771__model1_data4, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_ob10__0_053988__model1_data4 - 1.0*std::pow(-mob10__0_053988__model1_data4 + ob10__0_053988__model1_data4, 2)/std::pow(sigma_ob10__0_053988__model1_data4, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_ob10__0_057968__model1_data4 - 1.0*std::pow(-mob10__0_057968__model1_data4 + ob10__0_057968__model1_data4, 2)/std::pow(sigma_ob10__0_057968__model1_data4, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_ob10__0_060963__model1_data4 - 1.0*std::pow(-mob10__0_060963__model1_data4 + ob10__0_060963__model1_data4, 2)/std::pow(sigma_ob10__0_060963__model1_data4, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_ob10__0_066702__model1_data3 - 1.0*std::pow(-mob10__0_066702__model1_data3 + ob10__0_066702__model1_data3, 2)/std::pow(sigma_ob10__0_066702__model1_data3, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_ob10__0_068419__model1_data4 - 1.0*std::pow(-mob10__0_068419__model1_data4 + ob10__0_068419__model1_data4, 2)/std::pow(sigma_ob10__0_068419__model1_data4, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_ob10__0_092359__model1_data3 - 1.0*std::pow(-mob10__0_092359__model1_data3 + ob10__0_092359__model1_data3, 2)/std::pow(sigma_ob10__0_092359__model1_data3, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_ob10__0_10382__model1_data3 - 1.0*std::pow(-mob10__0_10382__model1_data3 + ob10__0_10382__model1_data3, 2)/std::pow(sigma_ob10__0_10382__model1_data3, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_ob10__0_110017__model1_data2 - 1.0*std::pow(-mob10__0_110017__model1_data2 + ob10__0_110017__model1_data2, 2)/std::pow(sigma_ob10__0_110017__model1_data2, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_ob10__0_127709__model1_data3 - 1.0*std::pow(-mob10__0_127709__model1_data3 + ob10__0_127709__model1_data3, 2)/std::pow(sigma_ob10__0_127709__model1_data3, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_ob10__0_138555__model1_data1 - 1.0*std::pow(-mob10__0_138555__model1_data1 + ob10__0_138555__model1_data1, 2)/std::pow(sigma_ob10__0_138555__model1_data1, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_ob10__0_147436__model1_data3 - 1.0*std::pow(-mob10__0_147436__model1_data3 + ob10__0_147436__model1_data3, 2)/std::pow(sigma_ob10__0_147436__model1_data3, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_ob10__0_167498__model1_data3 - 1.0*std::pow(-mob10__0_167498__model1_data3 + ob10__0_167498__model1_data3, 2)/std::pow(sigma_ob10__0_167498__model1_data3, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_ob10__0_194075__model1_data2 - 1.0*std::pow(-mob10__0_194075__model1_data2 + ob10__0_194075__model1_data2, 2)/std::pow(sigma_ob10__0_194075__model1_data2, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_ob10__0_223405__model1_data3 - 1.0*std::pow(-mob10__0_223405__model1_data3 + ob10__0_223405__model1_data3, 2)/std::pow(sigma_ob10__0_223405__model1_data3, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_ob10__0_234682__model1_data1 - 1.0*std::pow(-mob10__0_234682__model1_data1 + ob10__0_234682__model1_data1, 2)/std::pow(sigma_ob10__0_234682__model1_data1, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_ob10__0_240978__model1_data2 - 1.0*std::pow(-mob10__0_240978__model1_data2 + ob10__0_240978__model1_data2, 2)/std::pow(sigma_ob10__0_240978__model1_data2, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_ob10__0_249578__model1_data2 - 1.0*std::pow(-mob10__0_249578__model1_data2 + ob10__0_249578__model1_data2, 2)/std::pow(sigma_ob10__0_249578__model1_data2, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_ob10__0_267467__model1_data2 - 1.0*std::pow(-mob10__0_267467__model1_data2 + ob10__0_267467__model1_data2, 2)/std::pow(sigma_ob10__0_267467__model1_data2, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_ob10__0_29311__model1_data2 - 1.0*std::pow(-mob10__0_29311__model1_data2 + ob10__0_29311__model1_data2, 2)/std::pow(sigma_ob10__0_29311__model1_data2, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_ob10__0_298703__model1_data2 - 1.0*std::pow(-mob10__0_298703__model1_data2 + ob10__0_298703__model1_data2, 2)/std::pow(sigma_ob10__0_298703__model1_data2, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_ob10__0_31831__model1_data1 - 1.0*std::pow(-mob10__0_31831__model1_data1 + ob10__0_31831__model1_data1, 2)/std::pow(sigma_ob10__0_31831__model1_data1, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_ob10__0_320749__model1_data1 - 1.0*std::pow(-mob10__0_320749__model1_data1 + ob10__0_320749__model1_data1, 2)/std::pow(sigma_ob10__0_320749__model1_data1, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_ob10__0_320935__model1_data1 - 1.0*std::pow(-mob10__0_320935__model1_data1 + ob10__0_320935__model1_data1, 2)/std::pow(sigma_ob10__0_320935__model1_data1, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_ob10__0_406451__model1_data1 - 1.0*std::pow(-mob10__0_406451__model1_data1 + ob10__0_406451__model1_data1, 2)/std::pow(sigma_ob10__0_406451__model1_data1, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_ob10__0_406728__model1_data1 - 1.0*std::pow(-mob10__0_406728__model1_data1 + ob10__0_406728__model1_data1, 2)/std::pow(sigma_ob10__0_406728__model1_data1, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_obcar__0_054089__model1_data2 - 1.0*std::pow(-mobcar__0_054089__model1_data2 + obcar__0_054089__model1_data2, 2)/std::pow(sigma_obcar__0_054089__model1_data2, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_obcar__0_063285__model1_data3 - 1.0*std::pow(-mobcar__0_063285__model1_data3 + obcar__0_063285__model1_data3, 2)/std::pow(sigma_obcar__0_063285__model1_data3, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_obcar__0_064455__model1_data2 - 1.0*std::pow(-mobcar__0_064455__model1_data2 + obcar__0_064455__model1_data2, 2)/std::pow(sigma_obcar__0_064455__model1_data2, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_obcar__0_076702__model1_data3 - 1.0*std::pow(-mobcar__0_076702__model1_data3 + obcar__0_076702__model1_data3, 2)/std::pow(sigma_obcar__0_076702__model1_data3, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_obcar__0_076761__model1_data2 - 1.0*std::pow(-mobcar__0_076761__model1_data2 + obcar__0_076761__model1_data2, 2)/std::pow(sigma_obcar__0_076761__model1_data2, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_obcar__0_098916__model1_data3 - 1.0*std::pow(-mobcar__0_098916__model1_data3 + obcar__0_098916__model1_data3, 2)/std::pow(sigma_obcar__0_098916__model1_data3, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_obcar__0_111703__model1_data2 - 1.0*std::pow(-mobcar__0_111703__model1_data2 + obcar__0_111703__model1_data2, 2)/std::pow(sigma_obcar__0_111703__model1_data2, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_obcar__0_115071__model1_data3 - 1.0*std::pow(-mobcar__0_115071__model1_data3 + obcar__0_115071__model1_data3, 2)/std::pow(sigma_obcar__0_115071__model1_data3, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_obcar__0_124582__model1_data3 - 1.0*std::pow(-mobcar__0_124582__model1_data3 + obcar__0_124582__model1_data3, 2)/std::pow(sigma_obcar__0_124582__model1_data3, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_obcar__0_145321__model1_data2 - 1.0*std::pow(-mobcar__0_145321__model1_data2 + obcar__0_145321__model1_data2, 2)/std::pow(sigma_obcar__0_145321__model1_data2, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_obcar__0_146065__model1_data3 - 1.0*std::pow(-mobcar__0_146065__model1_data3 + obcar__0_146065__model1_data3, 2)/std::pow(sigma_obcar__0_146065__model1_data3, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_obcar__0_161073__model1_data3 - 1.0*std::pow(-mobcar__0_161073__model1_data3 + obcar__0_161073__model1_data3, 2)/std::pow(sigma_obcar__0_161073__model1_data3, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_obcar__0_179664__model1_data2 - 1.0*std::pow(-mobcar__0_179664__model1_data2 + obcar__0_179664__model1_data2, 2)/std::pow(sigma_obcar__0_179664__model1_data2, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_obcar__0_251006__model1_data2 - 1.0*std::pow(-mobcar__0_251006__model1_data2 + obcar__0_251006__model1_data2, 2)/std::pow(sigma_obcar__0_251006__model1_data2, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_obcry__0_137702__model1_data4 - 1.0*std::pow(-mobcry__0_137702__model1_data4 + obcry__0_137702__model1_data4, 2)/std::pow(sigma_obcry__0_137702__model1_data4, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_obcry__0_224241__model1_data4 - 1.0*std::pow(-mobcry__0_224241__model1_data4 + obcry__0_224241__model1_data4, 2)/std::pow(sigma_obcry__0_224241__model1_data4, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_obcry__0_241465__model1_data4 - 1.0*std::pow(-mobcry__0_241465__model1_data4 + obcry__0_241465__model1_data4, 2)/std::pow(sigma_obcry__0_241465__model1_data4, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_obcry__0_275261__model1_data4 - 1.0*std::pow(-mobcry__0_275261__model1_data4 + obcry__0_275261__model1_data4, 2)/std::pow(sigma_obcry__0_275261__model1_data4, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_obcry__0_338201__model1_data4 - 1.0*std::pow(-mobcry__0_338201__model1_data4 + obcry__0_338201__model1_data4, 2)/std::pow(sigma_obcry__0_338201__model1_data4, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_obcry__0_353484__model1_data4 - 1.0*std::pow(-mobcry__0_353484__model1_data4 + obcry__0_353484__model1_data4, 2)/std::pow(sigma_obcry__0_353484__model1_data4, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_obcry__0_43054__model1_data4 - 1.0*std::pow(-mobcry__0_43054__model1_data4 + obcry__0_43054__model1_data4, 2)/std::pow(sigma_obcry__0_43054__model1_data4, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_oohb10__0_023108__model1_data6 - 1.0*std::pow(-moohb10__0_023108__model1_data6 + oohb10__0_023108__model1_data6, 2)/std::pow(sigma_oohb10__0_023108__model1_data6, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_oohb10__0_030931__model1_data6 - 1.0*std::pow(-moohb10__0_030931__model1_data6 + oohb10__0_030931__model1_data6, 2)/std::pow(sigma_oohb10__0_030931__model1_data6, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_oohb10__0_032031__model1_data6 - 1.0*std::pow(-moohb10__0_032031__model1_data6 + oohb10__0_032031__model1_data6, 2)/std::pow(sigma_oohb10__0_032031__model1_data6, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_oohb10__0_041749__model1_data6 - 1.0*std::pow(-moohb10__0_041749__model1_data6 + oohb10__0_041749__model1_data6, 2)/std::pow(sigma_oohb10__0_041749__model1_data6, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_oohb10__0_064588__model1_data6 - 1.0*std::pow(-moohb10__0_064588__model1_data6 + oohb10__0_064588__model1_data6, 2)/std::pow(sigma_oohb10__0_064588__model1_data6, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_oohb10__0_076127__model1_data6 - 1.0*std::pow(-moohb10__0_076127__model1_data6 + oohb10__0_076127__model1_data6, 2)/std::pow(sigma_oohb10__0_076127__model1_data6, 3);
            break;
        case 55:
            dJydsigma[55] = 1.0/sigma_oohb10__0_090733__model1_data6 - 1.0*std::pow(-moohb10__0_090733__model1_data6 + oohb10__0_090733__model1_data6, 2)/std::pow(sigma_oohb10__0_090733__model1_data6, 3);
            break;
        case 56:
            dJydsigma[56] = 1.0/sigma_oohb10__0_117543__model1_data5 - 1.0*std::pow(-moohb10__0_117543__model1_data5 + oohb10__0_117543__model1_data5, 2)/std::pow(sigma_oohb10__0_117543__model1_data5, 3);
            break;
        case 57:
            dJydsigma[57] = 1.0/sigma_oohb10__0_13378__model1_data5 - 1.0*std::pow(-moohb10__0_13378__model1_data5 + oohb10__0_13378__model1_data5, 2)/std::pow(sigma_oohb10__0_13378__model1_data5, 3);
            break;
        case 58:
            dJydsigma[58] = 1.0/sigma_oohb10__0_146042__model1_data4 - 1.0*std::pow(-moohb10__0_146042__model1_data4 + oohb10__0_146042__model1_data4, 2)/std::pow(sigma_oohb10__0_146042__model1_data4, 3);
            break;
        case 59:
            dJydsigma[59] = 1.0/sigma_oohb10__0_184647__model1_data5 - 1.0*std::pow(-moohb10__0_184647__model1_data5 + oohb10__0_184647__model1_data5, 2)/std::pow(sigma_oohb10__0_184647__model1_data5, 3);
            break;
        case 60:
            dJydsigma[60] = 1.0/sigma_oohb10__0_185626__model1_data4 - 1.0*std::pow(-moohb10__0_185626__model1_data4 + oohb10__0_185626__model1_data4, 2)/std::pow(sigma_oohb10__0_185626__model1_data4, 3);
            break;
        case 61:
            dJydsigma[61] = 1.0/sigma_oohb10__0_191578__model1_data5 - 1.0*std::pow(-moohb10__0_191578__model1_data5 + oohb10__0_191578__model1_data5, 2)/std::pow(sigma_oohb10__0_191578__model1_data5, 3);
            break;
        case 62:
            dJydsigma[62] = 1.0/sigma_oohb10__0_205133__model1_data4 - 1.0*std::pow(-moohb10__0_205133__model1_data4 + oohb10__0_205133__model1_data4, 2)/std::pow(sigma_oohb10__0_205133__model1_data4, 3);
            break;
        case 63:
            dJydsigma[63] = 1.0/sigma_oohb10__0_208483__model1_data5 - 1.0*std::pow(-moohb10__0_208483__model1_data5 + oohb10__0_208483__model1_data5, 2)/std::pow(sigma_oohb10__0_208483__model1_data5, 3);
            break;
        case 64:
            dJydsigma[64] = 1.0/sigma_oohb10__0_249545__model1_data5 - 1.0*std::pow(-moohb10__0_249545__model1_data5 + oohb10__0_249545__model1_data5, 2)/std::pow(sigma_oohb10__0_249545__model1_data5, 3);
            break;
        case 65:
            dJydsigma[65] = 1.0/sigma_oohb10__0_259748__model1_data4 - 1.0*std::pow(-moohb10__0_259748__model1_data4 + oohb10__0_259748__model1_data4, 2)/std::pow(sigma_oohb10__0_259748__model1_data4, 3);
            break;
        case 66:
            dJydsigma[66] = 1.0/sigma_oohb10__0_270906__model1_data4 - 1.0*std::pow(-moohb10__0_270906__model1_data4 + oohb10__0_270906__model1_data4, 2)/std::pow(sigma_oohb10__0_270906__model1_data4, 3);
            break;
        case 67:
            dJydsigma[67] = 1.0/sigma_oohb10__0_319353__model1_data4 - 1.0*std::pow(-moohb10__0_319353__model1_data4 + oohb10__0_319353__model1_data4, 2)/std::pow(sigma_oohb10__0_319353__model1_data4, 3);
            break;
        case 68:
            dJydsigma[68] = 1.0/sigma_oohb10__0_365852__model1_data5 - 1.0*std::pow(-moohb10__0_365852__model1_data5 + oohb10__0_365852__model1_data5, 2)/std::pow(sigma_oohb10__0_365852__model1_data5, 3);
            break;
        case 69:
            dJydsigma[69] = 1.0/sigma_oohb10__0_387509__model1_data4 - 1.0*std::pow(-moohb10__0_387509__model1_data4 + oohb10__0_387509__model1_data4, 2)/std::pow(sigma_oohb10__0_387509__model1_data4, 3);
            break;
        case 70:
            dJydsigma[70] = 1.0/sigma_ozea__0_133086__model1_data6 - 1.0*std::pow(-mozea__0_133086__model1_data6 + ozea__0_133086__model1_data6, 2)/std::pow(sigma_ozea__0_133086__model1_data6, 3);
            break;
        case 71:
            dJydsigma[71] = 1.0/sigma_ozea__0_153161__model1_data6 - 1.0*std::pow(-mozea__0_153161__model1_data6 + ozea__0_153161__model1_data6, 2)/std::pow(sigma_ozea__0_153161__model1_data6, 3);
            break;
        case 72:
            dJydsigma[72] = 1.0/sigma_ozea__0_161417__model1_data6 - 1.0*std::pow(-mozea__0_161417__model1_data6 + ozea__0_161417__model1_data6, 2)/std::pow(sigma_ozea__0_161417__model1_data6, 3);
            break;
        case 73:
            dJydsigma[73] = 1.0/sigma_ozea__0_172674__model1_data6 - 1.0*std::pow(-mozea__0_172674__model1_data6 + ozea__0_172674__model1_data6, 2)/std::pow(sigma_ozea__0_172674__model1_data6, 3);
            break;
        case 74:
            dJydsigma[74] = 1.0/sigma_ozea__0_172799__model1_data6 - 1.0*std::pow(-mozea__0_172799__model1_data6 + ozea__0_172799__model1_data6, 2)/std::pow(sigma_ozea__0_172799__model1_data6, 3);
            break;
        case 75:
            dJydsigma[75] = 1.0/sigma_ozea__0_17774__model1_data6 - 1.0*std::pow(-mozea__0_17774__model1_data6 + ozea__0_17774__model1_data6, 2)/std::pow(sigma_ozea__0_17774__model1_data6, 3);
            break;
        case 76:
            dJydsigma[76] = 1.0/sigma_ozea__0_183243__model1_data6 - 1.0*std::pow(-mozea__0_183243__model1_data6 + ozea__0_183243__model1_data6, 2)/std::pow(sigma_ozea__0_183243__model1_data6, 3);
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
