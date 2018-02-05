/* Please do not change these definitions! */
#define nObs      399
#define nCov        2
#define nLoci       1
#define ARE_Ys_PRESORTED 0
#define TRANSFORMED 0
#define CENSORED    1
#define THRESH    log(300.)/log2
#define CALCL2B     0
#define CALCDVAR    0

/* NB: data have been reformatted from the original, to facilitate
 * analysis.  Because allele ids are stored as integers, rather than
 * strings, the following representation is used:
 * A1S->1, A2S->2, A3S->3, A24S->4
 *
 * Ordering of values in this & the next vector (y) must be preserved!
 */

int X[nObs][nCov]={
2,2,
2,2,
3,2,
3,2,
24,2,
2,2,
3,2,
1,2,
24,1,
3,24,
1,1,
3,1,
3,24,
24,3,
2,1,
2,1,
24,3,
24,3,
24,1,
24,2,
1,2,
3,1,
2,1,
3,3,
3,1,
24,2,
3,2,
3,2,
24,2,
2,1,
3,2,
3,1,
1,24,
2,1,
1,3,
2,1,
2,1,
3,2,
2,2,
24,2,
2,1,
3,3,
3,1,
1,1,
3,1,
3,3,
1,1,
1,1,
1,2,
3,1,
3,2,
3,1,
3,2,
24,24,
24,2,
1,1,
2,1,
2,2,
1,1,
24,2,
24,2,
3,2,
3,1,
3,2,
24,3,
3,2,
3,2,
3,2,
2,1,
2,1,
2,3,
24,1,
2,1,
1,3,
2,1,
3,2,
2,1,
2,2,
1,1,
1,1,
1,3,
24,3,
2,2,
1,2,
24,3,
1,1,
24,2,
24,3,
3,1,
1,1,
1,3,
1,3,
3,2,
1,3,
2,2,
3,2,
3,2,
1,2,
3,24,
1,3,
3,3,
3,2,
3,24,
3,3,
3,2,
2,2,
1,2,
1,2,
3,3,
1,1,
1,1,
3,3,
3,1,
1,1,
1,1,
3,1,
1,1,
24,1,
2,2,
24,3,
1,1,
3,3,
3,2,
2,1,
3,2,
3,2,
1,24,
3,3,
2,2,
3,2,
3,2,
2,1,
3,2,
2,1,
3,2,
3,3,
24,2,
2,1,
3,2,
3,3,
2,1,
3,2,
3,1,
3,2,
24,1,
3,2,
2,1,
3,3,
2,2,
3,2,
1,3,
3,2,
3,3,
1,1,
24,3,
2,2,
3,2,
2,1,
24,2,
24,3,
2,1,
3,2,
2,1,
1,3,
2,2,
3,3,
3,2,
3,3,
3,2,
1,2,
1,3,
2,2,
2,1,
2,1,
3,1,
1,24,
1,3,
3,2,
3,2,
24,2,
2,2,
24,2,
2,1,
3,3,
24,24,
3,3,
2,2,
3,2,
1,1,
24,1,
3,1,
24,2,
3,2,
2,1,
1,3,
2,1,
2,2,
24,2,
2,1,
24,2,
2,1,
24,24,
3,24,
3,1,
3,3,
24,3,
1,2,
1,1,
2,2,
2,1,
24,1,
3,2,
3,2,
24,1,
3,1,
3,1,
2,2,
3,3,
24,3,
24,2,
3,1,
3,2,
1,3,
3,24,
1,2,
1,3,
24,1,
3,2,
3,3,
2,1,
2,2,
1,1,
1,1,
2,2,
3,2,
2,2,
1,1,
1,2,
24,1,
3,3,
24,1,
3,1,
2,2,
3,2,
1,1,
3,3,
3,24,
24,24,
3,2,
2,1,
2,1,
3,1,
3,2,
3,3,
1,2,
3,1,
24,3,
2,1,
3,3,
3,1,
3,2,
3,2,
24,1,
24,3,
2,2,
2,1,
2,1,
3,24,
24,2,
3,1,
3,2,
3,2,
24,2,
2,2,
2,2,
2,1,
1,24,
24,3,
3,1,
1,1,
24,2,
1,2,
24,3,
24,3,
2,2,
3,2,
24,2,
3,3,
3,2,
3,1,
24,2,
1,1,
3,1,
24,3,
24,3,
1,3,
1,2,
24,2,
2,1,
1,1,
3,2,
3,2,
3,24,
24,2,
2,2,
2,1,
3,2,
24,24,
24,2,
24,2,
24,1,
2,2,
24,2,
24,2,
3,2,
2,2,
24,1,
3,1,
3,3,
24,3,
1,3,
3,1,
24,1,
3,1,
3,2,
2,24,
1,2,
3,2,
3,3,
2,1,
2,2,
3,24,
24,2,
2,2,
3,3,
2,2,
24,2,
24,2,
2,24,
3,2,
24,24,
3,2,
3,3,
2,2,
3,2,
3,1,
2,2,
24,1,
24,3,
2,2,
24,1,
3,2,
24,1,
3,2,
24,1,
3,24,
24,3,
2,1,
24,3,
2,24,
24,3,
3,1,
2,24,
3,24,
3,1,
24,24,
2,3,
24,1,
2,2,
24,3,
24,2,
3,2,
2,1,
1,2,
24,2,
2,2,
3,2,
24,2,
24,3,
3,2,
24,2,
3,2,
24,2,
3,3,
24,1,
1,2,
3,2,
24,2,
2,2,
24,24,
24,1,
3,3,
24,3,
3,24,
2,2,
24,3,
2,2,
3,2,
3,2
};

double y[nObs]={ 17930., 18735., 11529., 87328., 22207., 34008.,
11460., 532., 3898., 1076., 18380., 400., 400., 853., 87483., 13536.,
18880., 31146., 15200., 90390., 34780., 300., 9082., 20597., 1114.,
84442., 75667., 17311., 80827., 4810., 16270., 1117., 73376., 73777.,
300., 300., 19166., 3726., 38275., 11178., 10360., 15524., 1487.,
125423., 39334., 33999., 15948., 71021., 32536., 25220., 23474., 959.,
56599., 6157., 1786., 1058., 2385., 67428., 165471., 28002., 8047.,
150441., 74519., 5221., 50796., 30968., 736300., 43480., 154100.,
300., 5493., 23530., 23787., 21738., 2224., 3145., 79298., 7406.,
64391., 27900., 5490., 85160., 58266., 40493., 823542., 1894., 16802.,
63990., 772., 21299., 15350., 58817., 2145., 87681., 482., 89600.,
187323., 1208., 71789., 837., 44960., 524756., 24810., 168629., 3763.,
35350., 15281., 21108., 300., 39610., 37304., 4918., 33287., 13546.,
8095., 93248., 18700., 180556., 9741., 33941., 5785., 46074., 3034.,
60048., 45833., 20064., 8711., 1006., 519280., 2576., 4735., 6886.,
1118364., 2240., 11384., 20040., 34024., 10510., 36038., 25178.,
78619., 78426., 1487., 24377., 21047., 4694., 59679., 300., 119449.,
2511., 9817., 1050., 39283., 8537., 384246., 18351., 300., 12830.,
101523., 249300., 2120., 7651., 6330., 181514., 11060., 3546., 32050.,
20730., 9232., 1119., 1934., 12934., 4951., 18612., 105319., 7903.,
13263., 243901., 15657., 316823., 57375., 228401., 7539., 7755., 545.,
25480., 134681., 16681., 7433., 13950., 99370., 3429., 15363., 12177.,
8919., 17223., 515., 7647., 19190., 243750., 40098., 10046., 300.,
898., 2578., 73636., 12562., 6398., 300., 13615., 988., 99042.,
53881., 12810., 300., 2062., 12375., 9555., 56244., 51865., 28124.,
64240., 54053., 46750., 190000., 36706., 8511., 57632., 647., 23319.,
52840., 16606., 2186., 559243., 110021., 67256., 38350., 8262.,
44633., 908., 5580., 36726., 1309., 37788., 431., 53719., 2007.,
331927., 19373., 22614., 18590., 15056., 2530., 23266., 924., 957.,
12523., 57802., 8065., 31850., 8737., 13907., 15460., 1077., 13259.,
739., 10696., 4739., 2692., 703., 6835., 5215., 3527., 177481., 400.,
13349., 246943., 400., 40300., 18569., 1761., 4460., 59100., 14785.,
49889., 2552., 69923., 31660., 27599., 26507., 13250., 51131., 4950.,
11790., 13244., 16887., 7455., 6611., 4607., 6661., 40346., 6515.,
1529., 24826., 16338., 7043., 135473., 12950., 1874., 205914.,
109100., 20710., 1722., 400., 31912., 68640., 31121., 47640., 12510.,
1977., 126900., 4986., 86182., 19040., 34153., 91411., 7440., 8238.,
112700., 16420., 27820., 8403., 6377., 9748., 13460., 3434., 30208.,
12430., 3842., 1506806., 3066., 962., 31420., 29800., 28916., 11296.,
8990., 2706., 4137., 143611., 3845., 81658., 65455., 73035., 3098.,
10239., 7505., 104402., 593763., 35806., 4703., 300., 5557., 3083.,
300., 23928., 300., 32026., 718., 44551., 109556., 9496., 14024.,
17330., 18042., 14301., 625., 8764., 6648., 20856., 81425., 300.,
44564., 810., 6246., 35220., 948., 2741., 4564., 7898., 571., 27275.,
300., 18387., 6968., 78752., 4803., 156506., 71255.  };
