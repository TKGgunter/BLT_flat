#define muresolution_cxx
#include "BLT/BLTAnalysis/interface/muresolution.h"

const double muresolution::mtrk[12][13] = 
  {
    {  0.000000,  0.006387,  0.021692,  0.037723,  0.056343,  0.092962,  0.228220,  0.523078,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.005911,  0.021341,  0.038315,  0.057878,  0.092390,  0.202991,  0.485139,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.005746,  0.022573,  0.038588,  0.055072,  0.084012,  0.186516,  0.477003,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.005107,  0.018794,  0.034314,  0.051426,  0.080987,  0.189141,  0.490593,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.006285,  0.022758,  0.043502,  0.083725,  0.192718,  0.436864,  0.777533,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.013993,  0.034292,  0.068508,  0.180773,  0.480391,  0.855808,  0.999490,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.015649,  0.038022,  0.066657,  0.126886,  0.304814,  0.695729,  0.946563,  0.998691,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.005965,  0.016860,  0.036243,  0.057966,  0.094440,  0.183270,  0.349283,  0.582965,  0.861751,  0.984690,  0.999931,  1.000000},
    {  0.000000,  0.000953,  0.001922,  0.004746,  0.011688,  0.027296,  0.048902,  0.077316,  0.147987,  0.363561,  0.782071,  0.972489,  1.000000},
    {  0.000000,  0.000371,  0.000793,  0.001358,  0.002660,  0.008358,  0.031266,  0.104737,  0.275672,  0.605221,  0.923333,  0.999940,  1.000000},
    {  0.000000,  0.001940,  0.003245,  0.004432,  0.012021,  0.073583,  0.171702,  0.439224,  0.807839,  0.926306,  0.998857,  1.000000,  1.000000},
    {  0.000000,  0.014702,  0.071636,  0.103206,  0.157635,  0.273544,  0.476584,  0.926882,  0.971149,  0.999635,  0.999989,  1.000000,  1.000000}
  };
const double muresolution::dtrk[12][13] = 
  {
    {  0.000000,  0.015337,  0.042531,  0.066519,  0.093444,  0.152472,  0.309295,  0.607678,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.012074,  0.037115,  0.060476,  0.085501,  0.132165,  0.260040,  0.558338,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.012044,  0.038454,  0.060141,  0.081160,  0.120206,  0.240150,  0.550763,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.011261,  0.034859,  0.056024,  0.077894,  0.119067,  0.248101,  0.566697,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.012139,  0.037527,  0.066264,  0.116380,  0.240650,  0.494266,  0.811105,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.022272,  0.053426,  0.102620,  0.237061,  0.548374,  0.888598,  0.999871,  1.000000,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.020214,  0.047585,  0.085526,  0.175889,  0.391597,  0.782718,  0.970420,  0.999553,  1.000000,  1.000000,  1.000000,  1.000000},
    {  0.000000,  0.008410,  0.022270,  0.047309,  0.076764,  0.126225,  0.234795,  0.410378,  0.626922,  0.883026,  0.988090,  0.999955,  1.000000},
    {  0.000000,  0.002487,  0.005260,  0.010844,  0.023166,  0.045240,  0.071714,  0.111067,  0.206400,  0.438982,  0.812463,  0.974768,  1.000000},
    {  0.000000,  0.002628,  0.006726,  0.012207,  0.018655,  0.029107,  0.059974,  0.153381,  0.355815,  0.680486,  0.947136,  0.999957,  1.000000},
    {  0.000000,  0.007199,  0.015366,  0.027065,  0.047028,  0.092640,  0.213482,  0.503062,  0.858609,  0.951618,  0.999849,  1.000000,  1.000000},
    {  0.000000,  0.024096,  0.080354,  0.126425,  0.213504,  0.333944,  0.573026,  0.956701,  0.986892,  0.999981,  1.000000,  1.000000,  1.000000}
  };


const double muresolution::rmsA[12][12] = 
  {
    {  0.050726,  0.035646,  0.026464,  0.021660,  0.017669,  0.012909,  0.011375,  0.010404,  0.010404,  0.010404,  0.010404,  0.010404},
    {  0.052788,  0.035849,  0.026011,  0.021216,  0.017890,  0.013548,  0.012046,  0.011169,  0.011169,  0.011169,  0.011169,  0.011169},
    {  0.052627,  0.036126,  0.026231,  0.021649,  0.018645,  0.014414,  0.012791,  0.011943,  0.011943,  0.011943,  0.011943,  0.011943},
    {  0.052927,  0.036561,  0.026932,  0.022295,  0.019230,  0.015047,  0.013575,  0.012821,  0.012821,  0.012821,  0.012821,  0.012821},
    {  0.050928,  0.035651,  0.026984,  0.022483,  0.019761,  0.017719,  0.015704,  0.013821,  0.013821,  0.013821,  0.013821,  0.013821},
    {  0.051051,  0.034847,  0.025712,  0.021313,  0.019536,  0.018305,  0.016861,  0.016861,  0.016861,  0.016861,  0.016861,  0.016861},
    {  0.048635,  0.036267,  0.027832,  0.022982,  0.021264,  0.020433,  0.020130,  0.019816,  0.019816,  0.019816,  0.019816,  0.019816},
    {  0.057193,  0.043741,  0.032363,  0.029342,  0.025666,  0.023154,  0.021362,  0.020588,  0.019712,  0.018785,  0.018469,  0.018469},
    {  0.068717,  0.052634,  0.041876,  0.038327,  0.033866,  0.030875,  0.026821,  0.023181,  0.020975,  0.020349,  0.019712,  0.018880},
    {  0.083945,  0.076853,  0.062879,  0.052933,  0.041329,  0.035723,  0.031242,  0.028364,  0.025913,  0.025032,  0.024261,  0.024261},
    {  0.110282,  0.080547,  0.065664,  0.055678,  0.048897,  0.039287,  0.036640,  0.033058,  0.031977,  0.029760,  0.026956,  0.026956},
    {  0.109408,  0.087782,  0.072874,  0.064539,  0.056279,  0.045494,  0.041939,  0.040379,  0.038219,  0.038219,  0.038219,  0.038219}
  };
const double muresolution::rmsB[12][12] = 
  {
    {  0.000954,  0.000614,  0.000379,  0.000262,  0.000183,  0.000120,  0.000094,  0.000079,  0.000079,  0.000079,  0.000079,  0.000079},
    {  0.000992,  0.000620,  0.000335,  0.000219,  0.000176,  0.000114,  0.000088,  0.000072,  0.000072,  0.000072,  0.000072,  0.000072},
    {  0.000970,  0.000596,  0.000351,  0.000208,  0.000158,  0.000104,  0.000083,  0.000071,  0.000071,  0.000071,  0.000071,  0.000071},
    {  0.000915,  0.000592,  0.000327,  0.000215,  0.000151,  0.000107,  0.000081,  0.000070,  0.000070,  0.000070,  0.000070,  0.000070},
    {  0.000819,  0.000515,  0.000284,  0.000184,  0.000145,  0.000119,  0.000089,  0.000071,  0.000071,  0.000071,  0.000071,  0.000071},
    {  0.000812,  0.000399,  0.000202,  0.000118,  0.000097,  0.000090,  0.000083,  0.000083,  0.000083,  0.000083,  0.000083,  0.000083},
    {  0.000600,  0.000391,  0.000172,  0.000090,  0.000080,  0.000073,  0.000075,  0.000080,  0.000080,  0.000080,  0.000080,  0.000080},
    {  0.000770,  0.000448,  0.000213,  0.000137,  0.000107,  0.000113,  0.000091,  0.000090,  0.000092,  0.000086,  0.000085,  0.000085},
    {  0.000236,  0.000257,  0.000324,  0.000327,  0.000194,  0.000186,  0.000159,  0.000140,  0.000120,  0.000112,  0.000102,  0.000092},
    {  0.000756,  0.000760,  0.000577,  0.000156,  0.000282,  0.000307,  0.000256,  0.000210,  0.000169,  0.000154,  0.000140,  0.000140},
    {  0.001609,  0.000856,  0.000434,  0.000484,  0.000705,  0.000402,  0.000327,  0.000253,  0.000222,  0.000188,  0.000117,  0.000117},
    {  0.001966,  0.001466,  0.001063,  0.000993,  0.000835,  0.000559,  0.000467,  0.000363,  0.000307,  0.000307,  0.000307,  0.000307}
  };
const double muresolution::rmsC[12][12] = 
  {
    {  0.000000000,  0.000000818,  0.000001350,  0.000000818,  0.000000681,  0.000000496,  0.000000367,  0.000000290,  0.000000290,  0.000000290,  0.000000290,  0.000000290},
    {  0.000000000,  0.000001120,  0.000001624,  0.000001052,  0.000000453,  0.000000549,  0.000000444,  0.000000358,  0.000000358,  0.000000358,  0.000000358,  0.000000358},
    {  0.000000000,  0.000000672,  0.000000628,  0.000000545,  0.000000693,  0.000000539,  0.000000349,  0.000000364,  0.000000364,  0.000000364,  0.000000364,  0.000000364},
    {  0.000000000,  0.000000334,  0.000001170,  0.000000634,  0.000000466,  0.000000464,  0.000000350,  0.000000312,  0.000000312,  0.000000312,  0.000000312,  0.000000312},
    {  0.000000000,  0.000000883,  0.000000890,  0.000000809,  0.000000562,  0.000000457,  0.000000375,  0.000000284,  0.000000284,  0.000000284,  0.000000284,  0.000000284},
    {  0.000000000,  0.000001028,  0.000000631,  0.000000487,  0.000000369,  0.000000399,  0.000000269,  0.000000269,  0.000000269,  0.000000269,  0.000000269,  0.000000269},
    {  0.000001892,  0.000000973,  0.000000047,  0.000000553,  0.000000260,  0.000000172,  0.000000062,  0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000000000},
    {  0.000000000,  0.000001659,  0.000000899,  0.000000571,  0.000000591,  0.000000065,  0.000000212,  0.000000122,  0.000000158,  0.000000170,  0.000000222,  0.000000222},
    {  0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000001113,  0.000000927,  0.000000197,  0.000000432,  0.000000298,  0.000000270,  0.000000222,  0.000000000},
    {  0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000000378,  0.000000961,  0.000001055,  0.000000711,  0.000000559,  0.000000285,  0.000000285},
    {  0.000000000,  0.000000000,  0.000000000,  0.000000000,  0.000002156,  0.000001995,  0.000001756,  0.000001183,  0.000000532,  0.000001284,  0.000000000,  0.000000000},
    {  0.000000000,  0.000000000,  0.000000000,  0.000000987,  0.000002206,  0.000001655,  0.000001579,  0.000001324,  0.000000998,  0.000000998,  0.000000998,  0.000000998}
  };
const double muresolution::width[12][12] = 
  {
    {  0.773388,  0.929788,  0.891882,  0.873623,  0.789572,  0.850130,  0.913731,  0.943424,  0.943424,  0.943424,  0.943424,  0.943424},
    {  0.799082,  0.936089,  0.909605,  0.877136,  0.817363,  0.852899,  0.915383,  0.944236,  0.944236,  0.944236,  0.944236,  0.944236},
    {  0.790793,  0.938715,  0.905944,  0.880592,  0.843691,  0.876905,  0.925691,  0.946351,  0.946351,  0.946351,  0.946351,  0.946351},
    {  0.791327,  0.934117,  0.893207,  0.883378,  0.846871,  0.889917,  0.927519,  0.945721,  0.945721,  0.945721,  0.945721,  0.945721},
    {  0.808909,  0.911642,  0.892182,  0.902208,  0.917997,  0.927930,  0.939580,  0.949385,  0.949385,  0.949385,  0.949385,  0.949385},
    {  0.855318,  0.853253,  0.836964,  0.900230,  0.928380,  0.939471,  0.950064,  0.950064,  0.950064,  0.950064,  0.950064,  0.950064},
    {  0.841682,  0.834910,  0.868762,  0.901805,  0.925752,  0.939990,  0.945576,  0.948023,  0.948023,  0.948023,  0.948023,  0.948023},
    {  0.850933,  0.792941,  0.886669,  0.881651,  0.906239,  0.906871,  0.927462,  0.938668,  0.944354,  0.950895,  0.957422,  0.957422},
    {  0.866584,  0.808445,  0.842062,  0.880142,  0.904926,  0.899392,  0.895956,  0.896851,  0.928494,  0.942368,  0.946628,  0.950276},
    {  0.913471,  0.817522,  0.857430,  0.754740,  0.818824,  0.882116,  0.924401,  0.932482,  0.938191,  0.943957,  0.950292,  0.950292},
    {  0.891991,  0.916309,  0.912971,  0.784480,  0.853438,  0.931977,  0.933427,  0.937991,  0.940982,  0.951699,  1.018070,  1.018070},
    {  0.851430,  0.855598,  0.864440,  0.867223,  0.875031,  0.925754,  0.942935,  0.948304,  0.941438,  0.941438,  0.941438,  0.941438}
  };
const double muresolution::alpha[12][12] = 
  {
    {  1.227064,  2.054213,  2.052798,  1.896000,  1.511706,  1.710303,  1.988656,  2.191823,  2.191823,  2.191823,  2.191823,  2.191823},
    {  1.321400,  2.189076,  2.189183,  2.039016,  1.623196,  1.725041,  1.982250,  2.169758,  2.169758,  2.169758,  2.169758,  2.169758},
    {  1.332505,  2.145082,  2.139168,  2.105299,  1.836723,  1.875754,  2.095762,  2.191828,  2.191828,  2.191828,  2.191828,  2.191828},
    {  1.453816,  2.165156,  2.133490,  2.179228,  1.924729,  1.979393,  2.093592,  2.206372,  2.206372,  2.206372,  2.206372,  2.206372},
    {  1.521746,  2.085819,  2.061093,  2.138653,  2.109733,  2.071408,  2.128346,  2.240621,  2.240621,  2.240621,  2.240621,  2.240621},
    {  1.739504,  1.810347,  1.860896,  2.064940,  2.107142,  2.106725,  2.136929,  2.136929,  2.136929,  2.136929,  2.136929,  2.136929},
    {  1.679953,  1.636251,  1.968624,  2.125903,  2.176996,  2.245964,  2.230837,  2.274475,  2.274475,  2.274475,  2.274475,  2.274475},
    {  1.679512,  1.437708,  1.953868,  2.055867,  2.174302,  2.108087,  2.185773,  2.237150,  2.231316,  2.267243,  2.327168,  2.327168},
    {  1.786119,  1.545093,  2.050903,  2.192901,  2.324751,  2.257905,  2.174161,  2.189226,  2.214992,  2.271842,  2.296034,  2.366425},
    {  1.876599,  1.695806,  2.548160,  2.282824,  2.176284,  2.468697,  2.349571,  2.253769,  2.226360,  2.272395,  2.328231,  2.328231},
    {  2.206826,  2.353129,  2.138423,  2.303089,  2.426376,  2.316774,  2.241302,  2.279660,  2.281389,  2.357675,  2.544768,  2.544768},
    {  1.975418,  2.426682,  2.288963,  2.384072,  2.289668,  2.153149,  2.283300,  2.187157,  2.247811,  2.247811,  2.247811,  2.247811}
  };
const double muresolution::power[12][12] = 
  {
    {119.138597,  8.888120,  5.581096,  6.756484, 10.183593,  9.302869,  8.621751,  7.941881,  7.941881,  7.941881,  7.941881,  7.941881},
    { 47.707470,  6.782172,  4.795656,  5.111284,  8.903519,  9.157243,  9.057267,  8.805498,  8.805498,  8.805498,  8.805498,  8.805498},
    { 29.096709,  8.072847,  5.166171,  4.624730,  6.185009,  7.532592,  7.635325,  8.482485,  8.482485,  8.482485,  8.482485,  8.482485},
    { 11.677324,  6.932618,  4.876870,  4.127809,  5.270752,  6.613383,  7.867817,  7.984497,  7.984497,  7.984497,  7.984497,  7.984497},
    { 10.803078,  6.147882,  5.523705,  5.186823,  6.520073,  8.475438,  8.825668,  7.727975,  7.727975,  7.727975,  7.727975,  7.727975},
    {  8.509194,  6.987574,  5.577341,  5.948022,  7.624155,  9.460713, 11.411691, 11.411691, 11.411691, 11.411691, 11.411691, 11.411691},
    {  9.047397,  9.895516,  5.487944,  5.313541,  6.134606,  6.347844,  7.300583,  6.802426,  6.802426,  6.802426,  6.802426,  6.802426},
    {  9.404447, 13.939807,  6.690378,  5.214866,  5.002944,  5.772509,  6.105721,  6.372551,  7.065144,  7.597239,  7.204687,  7.204687},
    {  5.855467,  7.251417,  3.656077,  3.842618,  3.816364,  4.128978,  4.562428,  4.600460,  5.914993,  6.308324,  6.462586,  5.686825},
    {  3.349051,  3.964530,  1.500860,  1.979641,  3.075013,  2.814136,  4.442262,  5.715357,  6.556994,  6.454485,  6.164768,  6.164768},
    {  3.318704,  2.288174,  2.772171,  2.384296,  2.674520,  5.001897,  5.879100,  5.742804,  5.974316,  5.924937,  4.474561,  4.474561},
    {  4.644875,  2.554545,  3.238422,  2.945678,  3.438117,  6.663403,  6.156638,  8.597734,  6.293282,  6.293282,  6.293282,  6.293282}
  };
