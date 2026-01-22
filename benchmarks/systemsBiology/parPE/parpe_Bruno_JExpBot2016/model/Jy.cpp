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

void Jy_Bruno_JExpBot2016(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_043179__model1_data4, 2)) + 0.5*std::pow(-mob10__0_043179__model1_data4 + ob10__0_043179__model1_data4, 2)/std::pow(sigma_ob10__0_043179__model1_data4, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_046576__model1_data4, 2)) + 0.5*std::pow(-mob10__0_046576__model1_data4 + ob10__0_046576__model1_data4, 2)/std::pow(sigma_ob10__0_046576__model1_data4, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_048771__model1_data4, 2)) + 0.5*std::pow(-mob10__0_048771__model1_data4 + ob10__0_048771__model1_data4, 2)/std::pow(sigma_ob10__0_048771__model1_data4, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_053988__model1_data4, 2)) + 0.5*std::pow(-mob10__0_053988__model1_data4 + ob10__0_053988__model1_data4, 2)/std::pow(sigma_ob10__0_053988__model1_data4, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_057968__model1_data4, 2)) + 0.5*std::pow(-mob10__0_057968__model1_data4 + ob10__0_057968__model1_data4, 2)/std::pow(sigma_ob10__0_057968__model1_data4, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_060963__model1_data4, 2)) + 0.5*std::pow(-mob10__0_060963__model1_data4 + ob10__0_060963__model1_data4, 2)/std::pow(sigma_ob10__0_060963__model1_data4, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_066702__model1_data3, 2)) + 0.5*std::pow(-mob10__0_066702__model1_data3 + ob10__0_066702__model1_data3, 2)/std::pow(sigma_ob10__0_066702__model1_data3, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_068419__model1_data4, 2)) + 0.5*std::pow(-mob10__0_068419__model1_data4 + ob10__0_068419__model1_data4, 2)/std::pow(sigma_ob10__0_068419__model1_data4, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_092359__model1_data3, 2)) + 0.5*std::pow(-mob10__0_092359__model1_data3 + ob10__0_092359__model1_data3, 2)/std::pow(sigma_ob10__0_092359__model1_data3, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_10382__model1_data3, 2)) + 0.5*std::pow(-mob10__0_10382__model1_data3 + ob10__0_10382__model1_data3, 2)/std::pow(sigma_ob10__0_10382__model1_data3, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_110017__model1_data2, 2)) + 0.5*std::pow(-mob10__0_110017__model1_data2 + ob10__0_110017__model1_data2, 2)/std::pow(sigma_ob10__0_110017__model1_data2, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_127709__model1_data3, 2)) + 0.5*std::pow(-mob10__0_127709__model1_data3 + ob10__0_127709__model1_data3, 2)/std::pow(sigma_ob10__0_127709__model1_data3, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_138555__model1_data1, 2)) + 0.5*std::pow(-mob10__0_138555__model1_data1 + ob10__0_138555__model1_data1, 2)/std::pow(sigma_ob10__0_138555__model1_data1, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_147436__model1_data3, 2)) + 0.5*std::pow(-mob10__0_147436__model1_data3 + ob10__0_147436__model1_data3, 2)/std::pow(sigma_ob10__0_147436__model1_data3, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_167498__model1_data3, 2)) + 0.5*std::pow(-mob10__0_167498__model1_data3 + ob10__0_167498__model1_data3, 2)/std::pow(sigma_ob10__0_167498__model1_data3, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_194075__model1_data2, 2)) + 0.5*std::pow(-mob10__0_194075__model1_data2 + ob10__0_194075__model1_data2, 2)/std::pow(sigma_ob10__0_194075__model1_data2, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_223405__model1_data3, 2)) + 0.5*std::pow(-mob10__0_223405__model1_data3 + ob10__0_223405__model1_data3, 2)/std::pow(sigma_ob10__0_223405__model1_data3, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_234682__model1_data1, 2)) + 0.5*std::pow(-mob10__0_234682__model1_data1 + ob10__0_234682__model1_data1, 2)/std::pow(sigma_ob10__0_234682__model1_data1, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_240978__model1_data2, 2)) + 0.5*std::pow(-mob10__0_240978__model1_data2 + ob10__0_240978__model1_data2, 2)/std::pow(sigma_ob10__0_240978__model1_data2, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_249578__model1_data2, 2)) + 0.5*std::pow(-mob10__0_249578__model1_data2 + ob10__0_249578__model1_data2, 2)/std::pow(sigma_ob10__0_249578__model1_data2, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_267467__model1_data2, 2)) + 0.5*std::pow(-mob10__0_267467__model1_data2 + ob10__0_267467__model1_data2, 2)/std::pow(sigma_ob10__0_267467__model1_data2, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_29311__model1_data2, 2)) + 0.5*std::pow(-mob10__0_29311__model1_data2 + ob10__0_29311__model1_data2, 2)/std::pow(sigma_ob10__0_29311__model1_data2, 2);
            break;
        case 22:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_298703__model1_data2, 2)) + 0.5*std::pow(-mob10__0_298703__model1_data2 + ob10__0_298703__model1_data2, 2)/std::pow(sigma_ob10__0_298703__model1_data2, 2);
            break;
        case 23:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_31831__model1_data1, 2)) + 0.5*std::pow(-mob10__0_31831__model1_data1 + ob10__0_31831__model1_data1, 2)/std::pow(sigma_ob10__0_31831__model1_data1, 2);
            break;
        case 24:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_320749__model1_data1, 2)) + 0.5*std::pow(-mob10__0_320749__model1_data1 + ob10__0_320749__model1_data1, 2)/std::pow(sigma_ob10__0_320749__model1_data1, 2);
            break;
        case 25:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_320935__model1_data1, 2)) + 0.5*std::pow(-mob10__0_320935__model1_data1 + ob10__0_320935__model1_data1, 2)/std::pow(sigma_ob10__0_320935__model1_data1, 2);
            break;
        case 26:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_406451__model1_data1, 2)) + 0.5*std::pow(-mob10__0_406451__model1_data1 + ob10__0_406451__model1_data1, 2)/std::pow(sigma_ob10__0_406451__model1_data1, 2);
            break;
        case 27:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ob10__0_406728__model1_data1, 2)) + 0.5*std::pow(-mob10__0_406728__model1_data1 + ob10__0_406728__model1_data1, 2)/std::pow(sigma_ob10__0_406728__model1_data1, 2);
            break;
        case 28:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_054089__model1_data2, 2)) + 0.5*std::pow(-mobcar__0_054089__model1_data2 + obcar__0_054089__model1_data2, 2)/std::pow(sigma_obcar__0_054089__model1_data2, 2);
            break;
        case 29:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_063285__model1_data3, 2)) + 0.5*std::pow(-mobcar__0_063285__model1_data3 + obcar__0_063285__model1_data3, 2)/std::pow(sigma_obcar__0_063285__model1_data3, 2);
            break;
        case 30:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_064455__model1_data2, 2)) + 0.5*std::pow(-mobcar__0_064455__model1_data2 + obcar__0_064455__model1_data2, 2)/std::pow(sigma_obcar__0_064455__model1_data2, 2);
            break;
        case 31:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_076702__model1_data3, 2)) + 0.5*std::pow(-mobcar__0_076702__model1_data3 + obcar__0_076702__model1_data3, 2)/std::pow(sigma_obcar__0_076702__model1_data3, 2);
            break;
        case 32:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_076761__model1_data2, 2)) + 0.5*std::pow(-mobcar__0_076761__model1_data2 + obcar__0_076761__model1_data2, 2)/std::pow(sigma_obcar__0_076761__model1_data2, 2);
            break;
        case 33:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_098916__model1_data3, 2)) + 0.5*std::pow(-mobcar__0_098916__model1_data3 + obcar__0_098916__model1_data3, 2)/std::pow(sigma_obcar__0_098916__model1_data3, 2);
            break;
        case 34:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_111703__model1_data2, 2)) + 0.5*std::pow(-mobcar__0_111703__model1_data2 + obcar__0_111703__model1_data2, 2)/std::pow(sigma_obcar__0_111703__model1_data2, 2);
            break;
        case 35:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_115071__model1_data3, 2)) + 0.5*std::pow(-mobcar__0_115071__model1_data3 + obcar__0_115071__model1_data3, 2)/std::pow(sigma_obcar__0_115071__model1_data3, 2);
            break;
        case 36:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_124582__model1_data3, 2)) + 0.5*std::pow(-mobcar__0_124582__model1_data3 + obcar__0_124582__model1_data3, 2)/std::pow(sigma_obcar__0_124582__model1_data3, 2);
            break;
        case 37:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_145321__model1_data2, 2)) + 0.5*std::pow(-mobcar__0_145321__model1_data2 + obcar__0_145321__model1_data2, 2)/std::pow(sigma_obcar__0_145321__model1_data2, 2);
            break;
        case 38:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_146065__model1_data3, 2)) + 0.5*std::pow(-mobcar__0_146065__model1_data3 + obcar__0_146065__model1_data3, 2)/std::pow(sigma_obcar__0_146065__model1_data3, 2);
            break;
        case 39:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_161073__model1_data3, 2)) + 0.5*std::pow(-mobcar__0_161073__model1_data3 + obcar__0_161073__model1_data3, 2)/std::pow(sigma_obcar__0_161073__model1_data3, 2);
            break;
        case 40:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_179664__model1_data2, 2)) + 0.5*std::pow(-mobcar__0_179664__model1_data2 + obcar__0_179664__model1_data2, 2)/std::pow(sigma_obcar__0_179664__model1_data2, 2);
            break;
        case 41:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcar__0_251006__model1_data2, 2)) + 0.5*std::pow(-mobcar__0_251006__model1_data2 + obcar__0_251006__model1_data2, 2)/std::pow(sigma_obcar__0_251006__model1_data2, 2);
            break;
        case 42:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry__0_137702__model1_data4, 2)) + 0.5*std::pow(-mobcry__0_137702__model1_data4 + obcry__0_137702__model1_data4, 2)/std::pow(sigma_obcry__0_137702__model1_data4, 2);
            break;
        case 43:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry__0_224241__model1_data4, 2)) + 0.5*std::pow(-mobcry__0_224241__model1_data4 + obcry__0_224241__model1_data4, 2)/std::pow(sigma_obcry__0_224241__model1_data4, 2);
            break;
        case 44:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry__0_241465__model1_data4, 2)) + 0.5*std::pow(-mobcry__0_241465__model1_data4 + obcry__0_241465__model1_data4, 2)/std::pow(sigma_obcry__0_241465__model1_data4, 2);
            break;
        case 45:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry__0_275261__model1_data4, 2)) + 0.5*std::pow(-mobcry__0_275261__model1_data4 + obcry__0_275261__model1_data4, 2)/std::pow(sigma_obcry__0_275261__model1_data4, 2);
            break;
        case 46:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry__0_338201__model1_data4, 2)) + 0.5*std::pow(-mobcry__0_338201__model1_data4 + obcry__0_338201__model1_data4, 2)/std::pow(sigma_obcry__0_338201__model1_data4, 2);
            break;
        case 47:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry__0_353484__model1_data4, 2)) + 0.5*std::pow(-mobcry__0_353484__model1_data4 + obcry__0_353484__model1_data4, 2)/std::pow(sigma_obcry__0_353484__model1_data4, 2);
            break;
        case 48:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_obcry__0_43054__model1_data4, 2)) + 0.5*std::pow(-mobcry__0_43054__model1_data4 + obcry__0_43054__model1_data4, 2)/std::pow(sigma_obcry__0_43054__model1_data4, 2);
            break;
        case 49:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_023108__model1_data6, 2)) + 0.5*std::pow(-moohb10__0_023108__model1_data6 + oohb10__0_023108__model1_data6, 2)/std::pow(sigma_oohb10__0_023108__model1_data6, 2);
            break;
        case 50:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_030931__model1_data6, 2)) + 0.5*std::pow(-moohb10__0_030931__model1_data6 + oohb10__0_030931__model1_data6, 2)/std::pow(sigma_oohb10__0_030931__model1_data6, 2);
            break;
        case 51:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_032031__model1_data6, 2)) + 0.5*std::pow(-moohb10__0_032031__model1_data6 + oohb10__0_032031__model1_data6, 2)/std::pow(sigma_oohb10__0_032031__model1_data6, 2);
            break;
        case 52:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_041749__model1_data6, 2)) + 0.5*std::pow(-moohb10__0_041749__model1_data6 + oohb10__0_041749__model1_data6, 2)/std::pow(sigma_oohb10__0_041749__model1_data6, 2);
            break;
        case 53:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_064588__model1_data6, 2)) + 0.5*std::pow(-moohb10__0_064588__model1_data6 + oohb10__0_064588__model1_data6, 2)/std::pow(sigma_oohb10__0_064588__model1_data6, 2);
            break;
        case 54:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_076127__model1_data6, 2)) + 0.5*std::pow(-moohb10__0_076127__model1_data6 + oohb10__0_076127__model1_data6, 2)/std::pow(sigma_oohb10__0_076127__model1_data6, 2);
            break;
        case 55:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_090733__model1_data6, 2)) + 0.5*std::pow(-moohb10__0_090733__model1_data6 + oohb10__0_090733__model1_data6, 2)/std::pow(sigma_oohb10__0_090733__model1_data6, 2);
            break;
        case 56:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_117543__model1_data5, 2)) + 0.5*std::pow(-moohb10__0_117543__model1_data5 + oohb10__0_117543__model1_data5, 2)/std::pow(sigma_oohb10__0_117543__model1_data5, 2);
            break;
        case 57:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_13378__model1_data5, 2)) + 0.5*std::pow(-moohb10__0_13378__model1_data5 + oohb10__0_13378__model1_data5, 2)/std::pow(sigma_oohb10__0_13378__model1_data5, 2);
            break;
        case 58:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_146042__model1_data4, 2)) + 0.5*std::pow(-moohb10__0_146042__model1_data4 + oohb10__0_146042__model1_data4, 2)/std::pow(sigma_oohb10__0_146042__model1_data4, 2);
            break;
        case 59:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_184647__model1_data5, 2)) + 0.5*std::pow(-moohb10__0_184647__model1_data5 + oohb10__0_184647__model1_data5, 2)/std::pow(sigma_oohb10__0_184647__model1_data5, 2);
            break;
        case 60:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_185626__model1_data4, 2)) + 0.5*std::pow(-moohb10__0_185626__model1_data4 + oohb10__0_185626__model1_data4, 2)/std::pow(sigma_oohb10__0_185626__model1_data4, 2);
            break;
        case 61:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_191578__model1_data5, 2)) + 0.5*std::pow(-moohb10__0_191578__model1_data5 + oohb10__0_191578__model1_data5, 2)/std::pow(sigma_oohb10__0_191578__model1_data5, 2);
            break;
        case 62:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_205133__model1_data4, 2)) + 0.5*std::pow(-moohb10__0_205133__model1_data4 + oohb10__0_205133__model1_data4, 2)/std::pow(sigma_oohb10__0_205133__model1_data4, 2);
            break;
        case 63:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_208483__model1_data5, 2)) + 0.5*std::pow(-moohb10__0_208483__model1_data5 + oohb10__0_208483__model1_data5, 2)/std::pow(sigma_oohb10__0_208483__model1_data5, 2);
            break;
        case 64:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_249545__model1_data5, 2)) + 0.5*std::pow(-moohb10__0_249545__model1_data5 + oohb10__0_249545__model1_data5, 2)/std::pow(sigma_oohb10__0_249545__model1_data5, 2);
            break;
        case 65:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_259748__model1_data4, 2)) + 0.5*std::pow(-moohb10__0_259748__model1_data4 + oohb10__0_259748__model1_data4, 2)/std::pow(sigma_oohb10__0_259748__model1_data4, 2);
            break;
        case 66:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_270906__model1_data4, 2)) + 0.5*std::pow(-moohb10__0_270906__model1_data4 + oohb10__0_270906__model1_data4, 2)/std::pow(sigma_oohb10__0_270906__model1_data4, 2);
            break;
        case 67:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_319353__model1_data4, 2)) + 0.5*std::pow(-moohb10__0_319353__model1_data4 + oohb10__0_319353__model1_data4, 2)/std::pow(sigma_oohb10__0_319353__model1_data4, 2);
            break;
        case 68:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_365852__model1_data5, 2)) + 0.5*std::pow(-moohb10__0_365852__model1_data5 + oohb10__0_365852__model1_data5, 2)/std::pow(sigma_oohb10__0_365852__model1_data5, 2);
            break;
        case 69:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_oohb10__0_387509__model1_data4, 2)) + 0.5*std::pow(-moohb10__0_387509__model1_data4 + oohb10__0_387509__model1_data4, 2)/std::pow(sigma_oohb10__0_387509__model1_data4, 2);
            break;
        case 70:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea__0_133086__model1_data6, 2)) + 0.5*std::pow(-mozea__0_133086__model1_data6 + ozea__0_133086__model1_data6, 2)/std::pow(sigma_ozea__0_133086__model1_data6, 2);
            break;
        case 71:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea__0_153161__model1_data6, 2)) + 0.5*std::pow(-mozea__0_153161__model1_data6 + ozea__0_153161__model1_data6, 2)/std::pow(sigma_ozea__0_153161__model1_data6, 2);
            break;
        case 72:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea__0_161417__model1_data6, 2)) + 0.5*std::pow(-mozea__0_161417__model1_data6 + ozea__0_161417__model1_data6, 2)/std::pow(sigma_ozea__0_161417__model1_data6, 2);
            break;
        case 73:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea__0_172674__model1_data6, 2)) + 0.5*std::pow(-mozea__0_172674__model1_data6 + ozea__0_172674__model1_data6, 2)/std::pow(sigma_ozea__0_172674__model1_data6, 2);
            break;
        case 74:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea__0_172799__model1_data6, 2)) + 0.5*std::pow(-mozea__0_172799__model1_data6 + ozea__0_172799__model1_data6, 2)/std::pow(sigma_ozea__0_172799__model1_data6, 2);
            break;
        case 75:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea__0_17774__model1_data6, 2)) + 0.5*std::pow(-mozea__0_17774__model1_data6 + ozea__0_17774__model1_data6, 2)/std::pow(sigma_ozea__0_17774__model1_data6, 2);
            break;
        case 76:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ozea__0_183243__model1_data6, 2)) + 0.5*std::pow(-mozea__0_183243__model1_data6 + ozea__0_183243__model1_data6, 2)/std::pow(sigma_ozea__0_183243__model1_data6, 2);
            break;
    }
}

} // namespace model_Bruno_JExpBot2016
} // namespace amici
