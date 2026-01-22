#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <algorithm>

#include "p.h"
#include "y.h"
#include "sigmay.h"

namespace amici {
namespace model_Bruno_JExpBot2016 {

void sigmay_Bruno_JExpBot2016(realtype *sigmay, const realtype t, const realtype *p, const realtype *k, const realtype *y){
    sigma_ob10__0_043179__model1_data4 = noiseParameter1_ob10__0_043179__model1_data4;  // sigmay[0]
    sigma_ob10__0_046576__model1_data4 = noiseParameter1_ob10__0_046576__model1_data4;  // sigmay[1]
    sigma_ob10__0_048771__model1_data4 = noiseParameter1_ob10__0_048771__model1_data4;  // sigmay[2]
    sigma_ob10__0_053988__model1_data4 = noiseParameter1_ob10__0_053988__model1_data4;  // sigmay[3]
    sigma_ob10__0_057968__model1_data4 = noiseParameter1_ob10__0_057968__model1_data4;  // sigmay[4]
    sigma_ob10__0_060963__model1_data4 = noiseParameter1_ob10__0_060963__model1_data4;  // sigmay[5]
    sigma_ob10__0_066702__model1_data3 = noiseParameter1_ob10__0_066702__model1_data3;  // sigmay[6]
    sigma_ob10__0_068419__model1_data4 = noiseParameter1_ob10__0_068419__model1_data4;  // sigmay[7]
    sigma_ob10__0_092359__model1_data3 = noiseParameter1_ob10__0_092359__model1_data3;  // sigmay[8]
    sigma_ob10__0_10382__model1_data3 = noiseParameter1_ob10__0_10382__model1_data3;  // sigmay[9]
    sigma_ob10__0_110017__model1_data2 = noiseParameter1_ob10__0_110017__model1_data2;  // sigmay[10]
    sigma_ob10__0_127709__model1_data3 = noiseParameter1_ob10__0_127709__model1_data3;  // sigmay[11]
    sigma_ob10__0_138555__model1_data1 = noiseParameter1_ob10__0_138555__model1_data1;  // sigmay[12]
    sigma_ob10__0_147436__model1_data3 = noiseParameter1_ob10__0_147436__model1_data3;  // sigmay[13]
    sigma_ob10__0_167498__model1_data3 = noiseParameter1_ob10__0_167498__model1_data3;  // sigmay[14]
    sigma_ob10__0_194075__model1_data2 = noiseParameter1_ob10__0_194075__model1_data2;  // sigmay[15]
    sigma_ob10__0_223405__model1_data3 = noiseParameter1_ob10__0_223405__model1_data3;  // sigmay[16]
    sigma_ob10__0_234682__model1_data1 = noiseParameter1_ob10__0_234682__model1_data1;  // sigmay[17]
    sigma_ob10__0_240978__model1_data2 = noiseParameter1_ob10__0_240978__model1_data2;  // sigmay[18]
    sigma_ob10__0_249578__model1_data2 = noiseParameter1_ob10__0_249578__model1_data2;  // sigmay[19]
    sigma_ob10__0_267467__model1_data2 = noiseParameter1_ob10__0_267467__model1_data2;  // sigmay[20]
    sigma_ob10__0_29311__model1_data2 = noiseParameter1_ob10__0_29311__model1_data2;  // sigmay[21]
    sigma_ob10__0_298703__model1_data2 = noiseParameter1_ob10__0_298703__model1_data2;  // sigmay[22]
    sigma_ob10__0_31831__model1_data1 = noiseParameter1_ob10__0_31831__model1_data1;  // sigmay[23]
    sigma_ob10__0_320749__model1_data1 = noiseParameter1_ob10__0_320749__model1_data1;  // sigmay[24]
    sigma_ob10__0_320935__model1_data1 = noiseParameter1_ob10__0_320935__model1_data1;  // sigmay[25]
    sigma_ob10__0_406451__model1_data1 = noiseParameter1_ob10__0_406451__model1_data1;  // sigmay[26]
    sigma_ob10__0_406728__model1_data1 = noiseParameter1_ob10__0_406728__model1_data1;  // sigmay[27]
    sigma_obcar__0_054089__model1_data2 = noiseParameter1_obcar__0_054089__model1_data2;  // sigmay[28]
    sigma_obcar__0_063285__model1_data3 = noiseParameter1_obcar__0_063285__model1_data3;  // sigmay[29]
    sigma_obcar__0_064455__model1_data2 = noiseParameter1_obcar__0_064455__model1_data2;  // sigmay[30]
    sigma_obcar__0_076702__model1_data3 = noiseParameter1_obcar__0_076702__model1_data3;  // sigmay[31]
    sigma_obcar__0_076761__model1_data2 = noiseParameter1_obcar__0_076761__model1_data2;  // sigmay[32]
    sigma_obcar__0_098916__model1_data3 = noiseParameter1_obcar__0_098916__model1_data3;  // sigmay[33]
    sigma_obcar__0_111703__model1_data2 = noiseParameter1_obcar__0_111703__model1_data2;  // sigmay[34]
    sigma_obcar__0_115071__model1_data3 = noiseParameter1_obcar__0_115071__model1_data3;  // sigmay[35]
    sigma_obcar__0_124582__model1_data3 = noiseParameter1_obcar__0_124582__model1_data3;  // sigmay[36]
    sigma_obcar__0_145321__model1_data2 = noiseParameter1_obcar__0_145321__model1_data2;  // sigmay[37]
    sigma_obcar__0_146065__model1_data3 = noiseParameter1_obcar__0_146065__model1_data3;  // sigmay[38]
    sigma_obcar__0_161073__model1_data3 = noiseParameter1_obcar__0_161073__model1_data3;  // sigmay[39]
    sigma_obcar__0_179664__model1_data2 = noiseParameter1_obcar__0_179664__model1_data2;  // sigmay[40]
    sigma_obcar__0_251006__model1_data2 = noiseParameter1_obcar__0_251006__model1_data2;  // sigmay[41]
    sigma_obcry__0_137702__model1_data4 = noiseParameter1_obcry__0_137702__model1_data4;  // sigmay[42]
    sigma_obcry__0_224241__model1_data4 = noiseParameter1_obcry__0_224241__model1_data4;  // sigmay[43]
    sigma_obcry__0_241465__model1_data4 = noiseParameter1_obcry__0_241465__model1_data4;  // sigmay[44]
    sigma_obcry__0_275261__model1_data4 = noiseParameter1_obcry__0_275261__model1_data4;  // sigmay[45]
    sigma_obcry__0_338201__model1_data4 = noiseParameter1_obcry__0_338201__model1_data4;  // sigmay[46]
    sigma_obcry__0_353484__model1_data4 = noiseParameter1_obcry__0_353484__model1_data4;  // sigmay[47]
    sigma_obcry__0_43054__model1_data4 = noiseParameter1_obcry__0_43054__model1_data4;  // sigmay[48]
    sigma_oohb10__0_023108__model1_data6 = noiseParameter1_oohb10__0_023108__model1_data6;  // sigmay[49]
    sigma_oohb10__0_030931__model1_data6 = noiseParameter1_oohb10__0_030931__model1_data6;  // sigmay[50]
    sigma_oohb10__0_032031__model1_data6 = noiseParameter1_oohb10__0_032031__model1_data6;  // sigmay[51]
    sigma_oohb10__0_041749__model1_data6 = noiseParameter1_oohb10__0_041749__model1_data6;  // sigmay[52]
    sigma_oohb10__0_064588__model1_data6 = noiseParameter1_oohb10__0_064588__model1_data6;  // sigmay[53]
    sigma_oohb10__0_076127__model1_data6 = noiseParameter1_oohb10__0_076127__model1_data6;  // sigmay[54]
    sigma_oohb10__0_090733__model1_data6 = noiseParameter1_oohb10__0_090733__model1_data6;  // sigmay[55]
    sigma_oohb10__0_117543__model1_data5 = noiseParameter1_oohb10__0_117543__model1_data5;  // sigmay[56]
    sigma_oohb10__0_13378__model1_data5 = noiseParameter1_oohb10__0_13378__model1_data5;  // sigmay[57]
    sigma_oohb10__0_146042__model1_data4 = noiseParameter1_oohb10__0_146042__model1_data4;  // sigmay[58]
    sigma_oohb10__0_184647__model1_data5 = noiseParameter1_oohb10__0_184647__model1_data5;  // sigmay[59]
    sigma_oohb10__0_185626__model1_data4 = noiseParameter1_oohb10__0_185626__model1_data4;  // sigmay[60]
    sigma_oohb10__0_191578__model1_data5 = noiseParameter1_oohb10__0_191578__model1_data5;  // sigmay[61]
    sigma_oohb10__0_205133__model1_data4 = noiseParameter1_oohb10__0_205133__model1_data4;  // sigmay[62]
    sigma_oohb10__0_208483__model1_data5 = noiseParameter1_oohb10__0_208483__model1_data5;  // sigmay[63]
    sigma_oohb10__0_249545__model1_data5 = noiseParameter1_oohb10__0_249545__model1_data5;  // sigmay[64]
    sigma_oohb10__0_259748__model1_data4 = noiseParameter1_oohb10__0_259748__model1_data4;  // sigmay[65]
    sigma_oohb10__0_270906__model1_data4 = noiseParameter1_oohb10__0_270906__model1_data4;  // sigmay[66]
    sigma_oohb10__0_319353__model1_data4 = noiseParameter1_oohb10__0_319353__model1_data4;  // sigmay[67]
    sigma_oohb10__0_365852__model1_data5 = noiseParameter1_oohb10__0_365852__model1_data5;  // sigmay[68]
    sigma_oohb10__0_387509__model1_data4 = noiseParameter1_oohb10__0_387509__model1_data4;  // sigmay[69]
    sigma_ozea__0_133086__model1_data6 = noiseParameter1_ozea__0_133086__model1_data6;  // sigmay[70]
    sigma_ozea__0_153161__model1_data6 = noiseParameter1_ozea__0_153161__model1_data6;  // sigmay[71]
    sigma_ozea__0_161417__model1_data6 = noiseParameter1_ozea__0_161417__model1_data6;  // sigmay[72]
    sigma_ozea__0_172674__model1_data6 = noiseParameter1_ozea__0_172674__model1_data6;  // sigmay[73]
    sigma_ozea__0_172799__model1_data6 = noiseParameter1_ozea__0_172799__model1_data6;  // sigmay[74]
    sigma_ozea__0_17774__model1_data6 = noiseParameter1_ozea__0_17774__model1_data6;  // sigmay[75]
    sigma_ozea__0_183243__model1_data6 = noiseParameter1_ozea__0_183243__model1_data6;  // sigmay[76]
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
