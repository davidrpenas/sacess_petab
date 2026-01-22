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
namespace model_Chen_MSB2009 {

void dydp_Chen_MSB2009(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *tcl, const realtype *dtcldp, const realtype *spl, const realtype *sspl){
    switch(ip) {
        case 152:
            dydp[0] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__0_1__model1_data1, 2);
            break;
        case 153:
            dydp[1] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__0_1__model1_data2, 2);
            break;
        case 154:
            dydp[2] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__0_1__model1_data3, 2);
            break;
        case 155:
            dydp[3] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__0_1__model1_data4, 2);
            break;
        case 156:
            dydp[4] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__0_36__model1_data1, 2);
            break;
        case 157:
            dydp[5] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__0_69__model1_data1, 2);
            break;
        case 158:
            dydp[6] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__0_98__model1_data3, 2);
            break;
        case 159:
            dydp[7] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__1_5__model1_data1, 2);
            break;
        case 160:
            dydp[8] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__1_6__model1_data4, 2);
            break;
        case 161:
            dydp[9] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__1_8__model1_data2, 2);
            break;
        case 162:
            dydp[10] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__2_1__model1_data1, 2);
            break;
        case 163:
            dydp[11] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__2_3__model1_data1, 2);
            break;
        case 164:
            dydp[12] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__2_3__model1_data2, 2);
            break;
        case 165:
            dydp[13] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__2_6__model1_data2, 2);
            break;
        case 166:
            dydp[14] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__2_8__model1_data2, 2);
            break;
        case 167:
            dydp[15] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__2_9__model1_data2, 2);
            break;
        case 168:
            dydp[16] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__3_1__model1_data3, 2);
            break;
        case 169:
            dydp[17] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__3_2__model1_data2, 2);
            break;
        case 170:
            dydp[18] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__3_5__model1_data1, 2);
            break;
        case 171:
            dydp[19] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__3_6__model1_data4, 2);
            break;
        case 172:
            dydp[20] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__4_5__model1_data3, 2);
            break;
        case 173:
            dydp[21] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__4_6__model1_data2, 2);
            break;
        case 174:
            dydp[22] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__4_6__model1_data4, 2);
            break;
        case 175:
            dydp[23] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__4_7__model1_data2, 2);
            break;
        case 176:
            dydp[24] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__5_1__model1_data3, 2);
            break;
        case 177:
            dydp[25] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__5_6__model1_data3, 2);
            break;
        case 178:
            dydp[26] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__7_4__model1_data1, 2);
            break;
        case 179:
            dydp[27] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__7_7__model1_data4, 2);
            break;
        case 180:
            dydp[28] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__9_7__model1_data3, 2);
            break;
        case 181:
            dydp[29] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__11_1__model1_data3, 2);
            break;
        case 182:
            dydp[30] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__12_1__model1_data4, 2);
            break;
        case 183:
            dydp[31] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__13_1__model1_data4, 2);
            break;
        case 184:
            dydp[32] = (-c472 - c484 - c497 - c498)/std::pow(observableParameter1_AKT_PP__AKT_t__14_1__model1_data3, 2);
            break;
        case 185:
            dydp[33] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data1, 2);
            break;
        case 186:
            dydp[34] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data2, 2);
            break;
        case 187:
            dydp[35] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data3, 2);
            break;
        case 188:
            dydp[36] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__0_1__model1_data4, 2);
            break;
        case 189:
            dydp[37] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__2_4__model1_data1, 2);
            break;
        case 190:
            dydp[38] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__2_5__model1_data1, 2);
            break;
        case 191:
            dydp[39] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__4_7__model1_data1, 2);
            break;
        case 192:
            dydp[40] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__6_8__model1_data1, 2);
            break;
        case 193:
            dydp[41] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__7_4__model1_data1, 2);
            break;
        case 194:
            dydp[42] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__8_9__model1_data1, 2);
            break;
        case 195:
            dydp[43] = (-2*c100 - 2*c104 - c130 - c131 - c132 - 2*c136 - c148 - c149 - 2*c15 - c150 - c151 - c152 - c153 - c162 - c163 - c164 - c165 - c166 - c167 - 2*c17 - c171 - c172 - c173 - c174 - c175 - c176 - 2*c18 - c180 - c181 - c182 - c183 - c184 - c185 - c189 - 2*c19 - c190 - c191 - c192 - c193 - c194 - c195 - c196 - c197 - c198 - c199 - 2*c20 - c200 - c201 - c202 - c203 - c204 - c205 - c206 - c207 - c208 - c209 - 2*c21 - c210 - c211 - c212 - c213 - c214 - c215 - c216 - c217 - c218 - c219 - c220 - c221 - c222 - c223 - c224 - c225 - c226 - c227 - c228 - c229 - 2*c23 - c230 - c231 - c232 - c233 - c234 - c235 - c236 - c237 - c238 - c239 - c240 - c241 - c242 - c243 - c244 - c245 - c246 - c247 - c248 - c249 - 2*c25 - c250 - c251 - c252 - c253 - c254 - c255 - c256 - c257 - c258 - c259 - c260 - c261 - c262 - c263 - 2*c264 - c265 - c266 - c267 - 2*c27 - c281 - c282 - 2*c29 - 2*c32 - 2*c33 - 2*c34 - 2*c35 - 2*c36 - 2*c37 - c409 - c410 - c411 - c412 - 2*c415 - c416 - 2*c419 - 2*c420 - c427 - c428 - c429 - 2*c431 - 2*c432 - c433 - c434 - c435 - c437 - c438 - c440 - c445 - c446 - c447 - 2*c448 - c449 - c450 - c451 - c464 - c465 - c466 - 2*c483 - 2*c486 - 2*c489 - 2*c5 - 2*c63 - 2*c64 - 2*c65 - 2*c66 - 2*c67 - 2*c68 - 2*c7 - 2*c8 - 2*c88 - 2*c89 - 2*c90 - 2*c91 - 2*c92 - 2*c93 - 2*c94 - 2*c95 - 2*c96 - 2*c97 - 2*c98 - 2*c99)/std::pow(observableParameter1_ERB_B1_P_tot__EGFR_t__13_1__model1_data1, 2);
            break;
        case 196:
            dydp[44] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__0_1__model1_data1, 2);
            break;
        case 197:
            dydp[45] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__0_1__model1_data2, 2);
            break;
        case 198:
            dydp[46] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__0_1__model1_data3, 2);
            break;
        case 199:
            dydp[47] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__0_1__model1_data4, 2);
            break;
        case 200:
            dydp[48] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__0_46__model1_data1, 2);
            break;
        case 201:
            dydp[49] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__0_94__model1_data2, 2);
            break;
        case 202:
            dydp[50] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__0_94__model1_data3, 2);
            break;
        case 203:
            dydp[51] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_04__model1_data4, 2);
            break;
        case 204:
            dydp[52] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_1__model1_data4, 2);
            break;
        case 205:
            dydp[53] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_2__model1_data1, 2);
            break;
        case 206:
            dydp[54] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_2__model1_data2, 2);
            break;
        case 207:
            dydp[55] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_2__model1_data3, 2);
            break;
        case 208:
            dydp[56] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_5__model1_data1, 2);
            break;
        case 209:
            dydp[57] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_5__model1_data4, 2);
            break;
        case 210:
            dydp[58] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_7__model1_data4, 2);
            break;
        case 211:
            dydp[59] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_9__model1_data1, 2);
            break;
        case 212:
            dydp[60] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_9__model1_data2, 2);
            break;
        case 213:
            dydp[61] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__1_9__model1_data3, 2);
            break;
        case 214:
            dydp[62] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__2_2__model1_data4, 2);
            break;
        case 215:
            dydp[63] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__2_3__model1_data2, 2);
            break;
        case 216:
            dydp[64] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__2_3__model1_data3, 2);
            break;
        case 217:
            dydp[65] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__2_7__model1_data2, 2);
            break;
        case 218:
            dydp[66] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__2_7__model1_data3, 2);
            break;
        case 219:
            dydp[67] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__2_8__model1_data1, 2);
            break;
        case 220:
            dydp[68] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__2_9__model1_data1, 2);
            break;
        case 221:
            dydp[69] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__3_0__model1_data4, 2);
            break;
        case 222:
            dydp[70] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__3_1__model1_data1, 2);
            break;
        case 223:
            dydp[71] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__3_4__model1_data2, 2);
            break;
        case 224:
            dydp[72] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__3_4__model1_data3, 2);
            break;
        case 225:
            dydp[73] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__3_7__model1_data1, 2);
            break;
        case 226:
            dydp[74] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__4_0__model1_data2, 2);
            break;
        case 227:
            dydp[75] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__4_0__model1_data3, 2);
            break;
        case 228:
            dydp[76] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__4_3__model1_data2, 2);
            break;
        case 229:
            dydp[77] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__4_3__model1_data3, 2);
            break;
        case 230:
            dydp[78] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__4_3__model1_data4, 2);
            break;
        case 231:
            dydp[79] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__4_7__model1_data4, 2);
            break;
        case 232:
            dydp[80] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__6_8__model1_data2, 2);
            break;
        case 233:
            dydp[81] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__6_8__model1_data3, 2);
            break;
        case 234:
            dydp[82] = (-c101 - c102 - c431 - c432 - c433 - c434 - c435 - c437 - c438 - c440 - c474 - c475 - c477 - c478 - c480 - c481 - c59 - c61 - c83 - c84 - c95 - c96 - c97 - c98)/std::pow(observableParameter1_ERK_PP__ERK_t__8_2__model1_data4, 2);
            break;
    }
}

} // namespace model_Chen_MSB2009
} // namespace amici
