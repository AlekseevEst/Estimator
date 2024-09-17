#include <iostream>
#include <catch2/catch.hpp>
#include "imm.h"
#include "ukf.h"
// #include "ifilter.h"
// #include "initFilters.h"
using namespace Catch::Benchmark;


template <class M>
struct InitImmFilterTest
{
    struct InitIMMFilterCV
    {
        template <class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterCV, TypeArgs...> &filter,
                        const M &meas,
                        const M &measNoise)
        {

            filter.correctInfo.X << -208519.491248049,425.273419307277,156007.178510836,-537.892046379775,26427.9640461299,2170.91089446241;
            
            filter.correctInfo.P << 169092.379126566,377948.846346883,210714.997455627,478346.109034063,85241.4209925032,281546.206936399,
                                    377948.846346883,2083508.78573451,470359.505911582,2648243.34858566,283996.581634025,1267765.2905385,
                                    210714.997455627,470359.505911582,286390.39523831,637704.953645069,18827.7471983906,122744.872182853,
                                    478346.109034063,2648243.34858566,637704.953645069,3593758.45302526,120162.944428779,394818.984941201,
                                    85241.4209925032,283996.581634025,18827.7471983906,120162.944428779,567715.9487235,1486768.43942116,
                                    281546.206936399,1267765.2905385,122744.872182853,394818.984941201,1486768.43942116,7550611.65412092;

            filter.Q.resize(3,3);
            filter.Q << 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;

            filter.R = measNoise;

            filter.paramsSigmaPoints.alpha = 1e-3;
            filter.paramsSigmaPoints.beta = 2.0;
            filter.paramsSigmaPoints.kappa = 0.0;
            // инициализируем фильтр
        }
    };

    using TypeFilterCV = UnscentedKalmanFilter<M, InitIMMFilterCV, FuncConstVel<M>, FuncMeasSphVr<M>, FuncControlMatrix_XvXYvYZvZ<M>>;

    struct InitIMMFilterCTxy
    {
        template <class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterCTxy, TypeArgs...> &filter,
                        const M &meas,
                        const M &measNoise)
        {


            filter.correctInfo.X <<  -208421.133782935,820.027409144691,156147.324860707,26.4572927282898,-1.7203043514364,26352.902465481,1867.31710611064;

            filter.correctInfo.P << 211673.563164425,548315.263821949,271461.304500054,722546.969228855,-763.274698968697,52714.9687281008,149996.237113829,
                                    548315.263821949,2764910.51151843,713659.587738369,3626685.32280715,-2918.3772104135,153781.498493784,741086.650522075,
                                    271461.304500054,713659.587738369,372945.711693148,985385.321002178,-1107.6867053245,-27543.5797928967,-64805.7012818891,
                                    722546.969228855,3626685.32280715,985385.321002178,4989823.30013118,-4571.0472497858,-66173.8180647444,-358808.653586368,
                                    -763.274698968697,-2918.3772104135,-1107.6867053245,-4571.0472497858,31.7793901065311,587.422203345255,2362.68653623995,
                                    52714.9687281008,153781.498493784,-27543.5797928967,-66173.8180647444,587.422203345255,592552.847883242,1587220.38265364,
                                    149996.237113829,741086.650522075,-64805.7012818891,-358808.653586368,2362.68653623995,1587220.38265364,7956892.66736922;


            filter.Q.resize(4,4);
            filter.Q << 1.0, 0.0, 0.0, 0.0,
                        0.0, 1.0, 0.0, 0.0,
                        0.0, 0.0, 1.0, 0.0,
                        0.0, 0.0, 0.0, 1.0;

            filter.R = measNoise;

            filter.paramsSigmaPoints.alpha = 1e-3;
            filter.paramsSigmaPoints.beta = 2.0;
            filter.paramsSigmaPoints.kappa = 0.0;
        }
    };
    using TypeFilterCTxy = UnscentedKalmanFilter<M, InitIMMFilterCTxy, FuncConstTurnXY<M>, FuncMeasSphVr<M>, FuncControlMatrix_XvXYvYZvZW<M>>;

    struct InitIMMFilterCA
    {
        template <class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterCA, TypeArgs...> &filter,
                        const M &meas,
                        const M &measNoise)
        {


        filter.correctInfo.X << -208519.205711607,426.419072953101,-0.216158080583127,156007.58327308,-536.25552521745,0.162794641864398,26427.7474345561,2170.03956120518,0.0283001470407567;

        filter.correctInfo.P << 169205.378666363,378401.114232786,-82.6314244564966,210875.231903508,478995.287980985,62.0755871900824,85155.6164954445,281200.836199468,10.6017658311052,
                                378401.114232786,2085320.30500512,-306.344447931755,471001.819330289,2650847.93902894,243.762421082866,283652.970185939,1266382.20225638,41.7296951663936,
                                -82.6314244564966,-306.344447931755,102.046088781482,-122.559405794673,-505.146525611567,-0.0759778741905217,64.8781600641572,261.684243255828,-0.0115292952213456,
                                210875.231903508,471001.819330289,-122.559405794673,286616.975513115,638621.640477353,91.9914346855785,18706.2535003116,122255.834294644,15.4195278331535,
                                478995.287980985,2650847.93902894,-505.146525611567,638621.640477353,3597467.57213246,397.491635835657,119671.185539397,392839.589448836,63.4926967947048,
                                62.0755871900824,243.762421082866,-0.0759778741905217,91.9914346855785,397.491635835657,102.005714746229,-47.8737331356408,-193.087813011219,0.00888617797648842,
                                85155.6164954445,283652.970185939,64.8781600641572,18706.2535003116,119671.185539397,-47.8737331356408,567781.096276618,1487030.97291837,-6.63390135709578,
                                281200.836199468,1266382.20225638,261.684243255828,122255.834294644,392839.589448836,-193.087813011219,1487030.97291837,7551672.91432994,-8.49369947272032,
                                10.6017658311052,41.7296951663936,-0.0115292952213456,15.4195278331535,63.4926967947048,0.00888617797648842,-6.63390135709578,-8.49369947272032,101.948665214491;

        filter.Q.resize(3,3);
        filter.Q <<     1.0,   0.0,    0.0,  
                        0.0,   1.0,    0.0,
                        0.0,   0.0,    1.0;
                      


        filter.R = measNoise;

        filter.paramsSigmaPoints.alpha = 1e-3;
        filter.paramsSigmaPoints.beta = 2.0;
        filter.paramsSigmaPoints.kappa = 0.0;
        }
    };
    using TypeFilterCA = UnscentedKalmanFilter<M, InitIMMFilterCA, FuncConstAcceleration<M>, FuncMeasSphVr<M>, FuncControlMatrix_XvXaXYvYaYZvZaZ<M>>;


    struct InitIMMFilterBal
    {
        template <class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterBal, TypeArgs...> &filter,
                        const M &meas,
                        const M &measNoise)
        {

            filter.correctInfo.X << -208517.343435803,433.125412456341,156009.64388326,-526.164287411199,1.0029742634018,26426.0668930553,2161.01376342888;
            
            filter.correctInfo.P << 169800.355685803,380244.247321146,211726.335177038,481592.98069662,1.22637200156926,84591.8269890301,279346.174235904,
                                    380244.247321146,2089813.21184042,473674.829783065,2657293.41591475,5.27912896693937,281545.420435118,1259789.86255955,
                                    211726.335177038,473674.829783065,287825.727492025,642308.593094578,1.83409482140715,17917.0233248169,119780.286243308,
                                    481592.98069662,2657293.41591475,642308.593094578,3606084.43884313,8.6095700651371,116900.349205886,385084.529440626,
                                    1.22637200156926,5.27912896693937,1.83409482140715,8.6095700651371,1.00009878680555,-0.957248897088506,-4.32852212275974,
                                    84591.8269890301,281545.420435118,17917.0233248169,116900.349205886,-0.957248897088506,568214.24161534,1487799.84038955,
                                    279346.174235904,1259789.86255955,119780.286243308,385084.529440626,-4.32852212275974,1487799.84038955,7549342.06136858;

            filter.Q.resize(4,4);
            filter.Q << 900.,  0,    0,      0,
                          0,  900.,  0,      0,
                          0,   0,  0.0009,  0,
                          0,   0,   0,     900.;
                        

            filter.R = measNoise;

            filter.paramsSigmaPoints.alpha = 1e-3;
            filter.paramsSigmaPoints.beta = 2.0;
            filter.paramsSigmaPoints.kappa = 0.0;
            // инициализируем фильтр
        }
    };
    using TypeFilterBal = UnscentedKalmanFilter<M, InitIMMFilterBal, FuncBalreentry<M>, FuncMeasSphVr<M>, FuncControlMatrix_XvXYvYZvZBal<M>>;


    std::vector<std::shared_ptr<IFilter<M>>> filters; // <- если нужно проитерироваться по фильтрам

    std::shared_ptr<TypeFilterCV> filterCV; //<- для других моделей также
    std::shared_ptr<TypeFilterCTxy> filterCTxy;
    std::shared_ptr<TypeFilterCA> filterCA;
    std::shared_ptr<TypeFilterBal> filterBal;
    /*по аналогии для других моделей*/

    InitImmFilterTest() : filterCV{std::make_shared<TypeFilterCV>()},
                       filterCTxy{std::make_shared<TypeFilterCTxy>()},
                       filterCA{std::make_shared<TypeFilterCA>()},
                       filterBal{std::make_shared<TypeFilterBal>()}

    {
        filters.push_back(filterCV);
        filters.push_back(filterCTxy);
        filters.push_back(filterCA);
        filters.push_back(filterBal);
    }

    void operator()(IMM<M, InitImmFilterTest> & filterIMM,
                    const M &meas,
                    const M &measNoise)
                    
    {
        filterCV->Initialization(meas, measNoise);
        filterCTxy->Initialization(meas,measNoise);
        filterCA->Initialization(meas,measNoise);
        filterBal->Initialization(meas,measNoise);

        filterIMM.mu_i << 0.29442196992995,0.277097645716319,0.116243746451727,0.312236637902004;

        filterIMM.p_ij << 0.8,0.0666666666666667,0.0666666666666667,0.0666666666666667,
                          0.1,0.7,0.1,0.1,
                          0.1,0.1,0.7,0.1,
                          0.0333333333333333,0.0333333333333333,0.0333333333333333,0.9;


    }
};

TEST_CASE("imm_Predict")
{
    Eigen::MatrixXd meas(4,1);
    Eigen::MatrixXd measNoise(4,4);
    measNoise << 0.04,0,0,0,
                 0,0.04,0,0,
                 0,0,10000.,0,
                 0,0,0,625;


    IMM<Eigen::MatrixXd,InitImmFilterTest<Eigen::MatrixXd>> imm;
    imm.Initialization(meas, measNoise);
    double dt = 0.2;

    // PRINTM(imm.mu_i);
    // PRINTM(imm.p_ij);
    auto pred = imm.predict(dt);


    Eigen::MatrixXd expectedMixModelProbabilities (1,4);
    expectedMixModelProbabilities << 0.292530533026345,0.267220135378333,0.121988294193237,0.318261037402084;

    Eigen::MatrixXd expectedMixStateConstvel (6,1);
    expectedMixStateConstvel << -208518.372806711,429.755588134847,156008.766956284,-531.479911030052,26427.1082961814,2167.43056930405;


    Eigen::MatrixXd expectedMixStateCovarianceConstvel(6,6);
    expectedMixStateCovarianceConstvel << 169680.952882304,380300.393883794,211554.461578461,481714.821788146,84790.9692936867,279728.09005694,
                                          380300.393883794,2092895.13895725,473716.610325867,2661710.69228655,282193.738478825,1260490.22009067,
                                          211554.461578461,473716.610325867,287586.435113333,642501.509826996,18185.7773036283,120154.638937808,
                                          481714.821788146,2661710.69228655,642501.509826996,3612975.6686242,117585.651265279,384425.219724197,
                                          84790.9692936867,282193.738478825,18185.7773036283,117585.651265279,568059.889229878,1488150.99138618,
                                          279728.09005694,1260490.22009067,120154.638937808,384425.219724197,1488150.99138618,7556156.47241516;



    Eigen::MatrixXd expectedMixStateConstturn (7,1);
    expectedMixStateConstturn << -208432.418141584,774.723485150105,156131.234787061,-38.300761603085,-1.52203447094937,26361.509290865,1902.08626450005;

    Eigen::MatrixXd expectedMixStateCovarianceConstturn (7,7);
    expectedMixStateCovarianceConstturn <<  207761.274838047,532665.545017444,265879.294511013,700098.518578362,-692.480341298106,55702.5156057859,162090.649394978,
                                            532665.545017444,2702322.37316198,691302.011005603,3536710.9710522,-2650.9802212762,165736.07691672,789490.724935769,
                                            265879.294511013,691302.011005603,364993.308326047,953435.569101449,-1004.51236602898,-23284.2114386194,-47559.7680971341,
                                            700098.518578362,3536710.9710522,953435.569101449,4861492.71416074,-4142.78528047196,-49057.5575639907,-289497.2130678,
                                            -692.480341298106,-2650.9802212762,-1004.51236602898,-4142.78528047196,39.9437794968103,532.820038941398,2143.30017268391,
                                            55702.5156057859,165736.07691672,-23284.2114386194,-49057.5575639907,532.820038941398,590270.972870406,1577970.09897705,
                                            162090.649394978,789490.724935769,-47559.7680971341,-289497.2130678,2143.30017268391,1577970.09897705,7919365.7234292;

    Eigen::MatrixXd expectedMixStateConstacc (9,1);
    expectedMixStateConstacc<< -208518.273909398,430.153177951378,-0.207670070621582,156008.90666484,-530.913408538829,0.156402086295329,26427.0343762553,2167.1389508029,0.0271888680668359;
    
    Eigen::MatrixXd expectedMixStateCovarianceConstacc (9,9);
    expectedMixStateCovarianceConstacc << 169697.384817039,380366.771473945,-79.1931789526433,211577.017925024,481811.126956347,59.4922917898647,84779.0278649247,279680.811545591,10.1601252388846,
                                          380366.771473945,2093165.86822397,-293.539573540651,473808.216474288,2662104.32654313,233.606432991104,282145.768329112,1260300.0322586,39.9895439234118,
                                          -79.1931789526433,-293.539573540651,101.967506482599,-117.471962691829,-484.201234012892,-0.0743219475623182,62.1824692832777,250.806159988814,-0.011307346802426,
                                          211577.017925024,473808.216474288,-117.471962691829,287616.951971369,642631.408310808,88.1721702026633,18169.5070139186,120090.124903953,14.7780583926766,
                                          481811.126956347,2662104.32654313,-484.201234012892,642631.408310808,3613530.22975796,381.047573397316,117516.571453776,384150.354138225,60.8542445992885,
                                          59.4922917898647,233.606432991104,-0.0743219475623182,88.1721702026633,381.047573397316,101.927954942384,-45.8823225235012,-185.052055667945,0.00871104540648178,
                                          84779.0278649247,282145.768329112,62.1824692832777,18169.5070139186,117516.571453776,-45.8823225235012,568068.659778766,1488186.82517293,-6.3540166959762,
                                          279680.811545591,1260300.0322586,250.806159988814,120090.124903953,384150.354138225,-185.052055667945,1488186.82517293,7556307.84381271,-8.08130792975397,
                                          10.1601252388846,39.9895439234118,-0.011307346802426,14.7780583926766,60.8542445992885,0.00871104540648178,-6.3540166959762,-8.08130792975397,101.87217601447;

    
    
    
    Eigen::MatrixXd expectedMixStateBalreentry (7,1);
    expectedMixStateBalreentry << -208516.628014704,436.027631266079,156010.687175684,-522.035438878444,1.00285459365647,26425.5313180067,2158.93711889139;

    Eigen::MatrixXd expectedMixStateCovarianceBalreentry(7,7);
    expectedMixStateCovarianceBalreentry << 170193.607144776,381840.714815802,212287.560488261,483881.614064807,1.17498658144134,84296.0489885349,278134.278835395,
                                            381840.714815802,2096320.54269455,475953.884743078,2666628.59212603,5.05843810076393,280354.333088001,1254900.5866876,
                                            212287.560488261,475953.884743078,288626.200844817,645570.539071036,1.75732168163349,17494.4783427015,118044.615229228,
                                            483881.614064807,2666628.59212603,645570.539071036,3619430.45596574,8.25137708377678,115186.625592482,378023.110796311,
                                            1.17498658144134,5.05843810076393,1.75732168163349,8.25137708377678,1.00009515371838,-0.917205055749871,-4.14843568429346,
                                            84296.0489885349,280354.333088001,17494.4783427015,115186.625592482,-0.917205055749871,568439.841306354,1488750.70911053,
                                            278134.278835395,1254900.5866876,118044.615229228,378023.110796311,-4.14843568429346,1488750.70911053,7553402.37602857;

   
   
    CHECK((imm.mu_i - expectedMixModelProbabilities).norm() == Approx(0.0).margin(1e-9));
    bool condition = (imm.mu_i - expectedMixModelProbabilities).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.mu_i);
        PRINTM(expectedMixModelProbabilities);
    }


    CHECK((imm.stateMixed[0] - expectedMixStateConstvel).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.stateMixed[0] - expectedMixStateConstvel).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.stateMixed[0]);
        PRINTM(expectedMixStateConstvel);
    }

    CHECK((imm.covarianceMixed[0] - expectedMixStateCovarianceConstvel).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.covarianceMixed[0] - expectedMixStateCovarianceConstvel).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.covarianceMixed[0]);
        PRINTM(expectedMixStateCovarianceConstvel);
    }



    CHECK((imm.stateMixed[1] - expectedMixStateConstturn).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.stateMixed[1] - expectedMixStateConstturn).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.stateMixed[1]);
        PRINTM(expectedMixStateConstturn);
    }


    CHECK((imm.covarianceMixed[1] - expectedMixStateCovarianceConstturn).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.covarianceMixed[1] - expectedMixStateCovarianceConstturn).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.covarianceMixed[1]);
        PRINTM(expectedMixStateCovarianceConstturn);
    }

    CHECK((imm.stateMixed[2] - expectedMixStateConstacc).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.stateMixed[2] - expectedMixStateConstacc).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.stateMixed[2]);
        PRINTM(expectedMixStateConstacc);
    }

    CHECK((imm.covarianceMixed[2] - expectedMixStateCovarianceConstacc).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.covarianceMixed[2] - expectedMixStateCovarianceConstacc).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.covarianceMixed[2]);
        PRINTM(expectedMixStateCovarianceConstacc);
    }


    CHECK((imm.stateMixed[3] - expectedMixStateBalreentry).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.stateMixed[3] - expectedMixStateBalreentry).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.stateMixed[3]);
        PRINTM(expectedMixStateBalreentry);
    }

    CHECK((imm.covarianceMixed[3] - expectedMixStateCovarianceBalreentry).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.covarianceMixed[3] - expectedMixStateCovarianceBalreentry).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.covarianceMixed[3]);
        PRINTM(expectedMixStateCovarianceBalreentry);
    }



    Eigen::MatrixXd expectedPredictStateConstacc (9,1);
    expectedPredictStateConstacc << -208432.247394822,430.111643859071,-0.207670070577889,155902.727097752,-530.882128003242,0.156402086257365,26860.462706581,2167.14438850287,0.0271888680626418;
    
    Eigen::MatrixXd expectedPredictStateCovarianceConstacc (9,9);
    expectedPredictStateCovarianceConstacc << 405565.253286002,798966.905982783,-135.841743522719,509182.897675259,1014243.55066401,106.212091949935,247558.955162722,531749.465685699,18.1578078767681,
                                                798966.905982783,2093052.5710948,-272.946072244149,1006190.89117701,2662054.20461005,233.591568601587,534229.043267595,1260358.19094723,39.9872824540518,
                                                -135.841743522719,-272.946072244149,102.967506482599,-214.31369591833,-484.216098402559,-0.0743219475622864,112.34347512999,250.803898518972,-0.0113073468024275,
                                                509182.897675259,1006190.89117701,-214.31369591833,689217.340926803,1365378.36326959,166.440243977464,81055.7415758174,196921.884503362,26.949081533538,
                                                1014243.55066401,2662054.20461005,-484.216098402559,1365378.36326959,3613686.76590554,401.633164385832,194331.280855692,384125.514924652,60.8559868083708,
                                                106.212091949935,233.591568601587,-0.0743219475622864,166.440243977464,401.633164385832,102.927954942384,-82.8925594337691,-185.050313458125,0.00871104540648047,
                                                247558.955162722,534229.043267595,112.34347512999,81055.7415758174,194331.280855692,-82.8925594337691,1465595.42593778,2999447.04974204,-5.91283476513522,
                                                531749.465685699,1260358.19094723,250.803898518972,196921.884503362,384125.514924652,-185.050313458125,2999447.04974204,7556308.72617666,12.4931272732507,
                                                18.1578078767681,39.9872824540518,-0.0113073468024275,26.949081533538,60.8559868083708,0.00871104540648047,-5.91283476513522,12.4931272732507,102.87217601447;
    
    Eigen::MatrixXd expectedPredictStateConstvel (6,1);
    expectedPredictStateConstvel << -208432.421632005,429.755588024374,155902.471001581,-531.479910877958,26860.5944044653,2167.43056909479;
    
    Eigen::MatrixXd expectedPredictStateCovarianceConstvel(6,6);
    expectedPredictStateCovarianceConstvel << 405516.916395118,798879.425675349,509109.175695506,1014056.96024559,247594.943804385,531826.134075379,
                                                798879.425675349,2092895.17895725,1006058.74878168,2661710.69228654,534291.782496918,1260490.22009064,
                                                509109.175695506,1006058.74878168,689106.066193021,1365096.64754997,81110.8441354848,197039.682882747,
                                                1014056.96024559,2661710.69228654,1365096.64754997,3612975.70862418,194470.69521003,384425.219724165,
                                                247594.943804385,534291.782496918,81110.8441354848,194470.69521003,1465566.54508098,2999382.289869,
                                                531826.134075379,1260490.22009064,197039.682882747,384425.219724165,2999382.289869,7556156.51241509;
    
    Eigen::MatrixXd expectedPredictStateConstturn(7,1);
    expectedPredictStateConstturn << -208276.064289581,788.732267051918,156122.233256262,-51.7362758701388,-1.52203447113459,26741.9265430011,1902.08626481257;
    
    Eigen::MatrixXd expectedPredictStateCovarianceConstturn(7,7);
    expectedPredictStateCovarianceConstturn << 530378.117484034,1082242.67990564,685721.23884069,1400933.18302126,-1224.30262363314,152804.135624665,319864.952415688,
                                                1082242.67990564,2738970.11238496,1406427.5401029,3540293.36412852,-2667.03881834235,323203.762130261,788258.857450161,
                                                685721.23884069,1406427.5401029,938349.786573045,1910352.45246143,-1820.85707843142,-54099.4425414939,-105298.88003302,
                                                1400933.18302126,3540293.36412852,1910352.45246143,4801733.25050641,-4020.65260654048,-106075.513434436,-287893.082837211,
                                                -1224.30262363314,-2667.03881834235,-1820.85707843142,-4020.65260654048,39.9837794968104,961.480073478676,2143.30017268385,
                                                152804.135624665,323203.762130261,-54099.4425414939,-106075.513434436,961.480073478676,1538233.64179937,3161843.24766369,
                                                319864.952415688,788258.857450161,-105298.88003302,-287893.082837211,2143.30017268385,3161843.24766369,7919365.76342927;
    
    Eigen::MatrixXd expectedPredictStateBalreentry(7,1);
    expectedPredictStateBalreentry << -208429.451705492,435.735661658015,155906.301423818,-521.821925680537,1.00285459341874,26856.9129827778,2154.88162140126;
    
    Eigen::MatrixXd expectedPredictStateCovarianceBalreentry(7,7);
    expectedPredictStateCovarianceBalreentry << 406724.275024456,800738.148962902,510854.335395877,1016805.7415674,2.17193086138893,246124.779411838,528726.919927095,
                                                800738.148962902,2094831.03008298,1008854.28732918,2664906.37803617,4.91100473131214,530913.942372049,1253322.29123371,
                                                510854.335395877,1008854.28732918,691540.88986038,1368883.32167066,3.42476800326222,79267.3323477087,193606.180713761,
                                                1016805.7415674,2664906.37803617,1368883.32167066,3617038.76869894,8.42308615315357,190918.283789354,378294.768535776,
                                                2.17193086138893,4.91100473131214,3.42476800326222,8.42308615315357,1.00013115371836,-1.81891141351345,-4.86862809193442,
                                                246124.779411838,530913.942372049,79267.3323477087,190918.283789354,-1.81891141351345,1465706.08147528,2997108.4984913,
                                                528726.919927095,1253322.29123371,193606.180713761,378294.768535776,-4.86862809193442,2997108.4984913,7543962.53078896;



    CHECK((imm.initializator.filterCV->getPredictInfo().Xe - expectedPredictStateConstvel).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterCV->getPredictInfo().Xe - expectedPredictStateConstvel).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterCV->getPredictInfo().Xe);
        PRINTM(expectedPredictStateConstvel);
    }
    CHECK((imm.initializator.filterCV->getPredictInfo().Pe - expectedPredictStateCovarianceConstvel).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterCV->getPredictInfo().Pe - expectedPredictStateCovarianceConstvel).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterCV->getPredictInfo().Pe);
        PRINTM(expectedPredictStateCovarianceConstvel);
    }

    CHECK((imm.initializator.filterCTxy->getPredictInfo().Xe - expectedPredictStateConstturn).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterCTxy->getPredictInfo().Xe - expectedPredictStateConstturn).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterCTxy->getPredictInfo().Xe);
        PRINTM(expectedPredictStateConstturn);
    }
    CHECK((imm.initializator.filterCTxy->getPredictInfo().Pe - expectedPredictStateCovarianceConstturn).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterCTxy->getPredictInfo().Pe - expectedPredictStateCovarianceConstturn).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterCTxy->getPredictInfo().Pe);
        PRINTM(expectedPredictStateCovarianceConstturn);
    }

    CHECK((imm.initializator.filterCA->getPredictInfo().Xe - expectedPredictStateConstacc).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterCA->getPredictInfo().Xe - expectedPredictStateConstacc).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterCA->getPredictInfo().Xe);
        PRINTM(expectedPredictStateConstacc);
    }

    CHECK((imm.initializator.filterCA->getPredictInfo().Pe - expectedPredictStateCovarianceConstacc).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterCA->getPredictInfo().Pe - expectedPredictStateCovarianceConstacc).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterCA->getPredictInfo().Pe);
        PRINTM(expectedPredictStateCovarianceConstacc);
    }

    CHECK((imm.initializator.filterBal->getPredictInfo().Xe - expectedPredictStateBalreentry).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterBal->getPredictInfo().Xe - expectedPredictStateBalreentry).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterBal->getPredictInfo().Xe);
        PRINTM(expectedPredictStateBalreentry);
    }
    CHECK((imm.initializator.filterBal->getPredictInfo().Pe - expectedPredictStateCovarianceBalreentry).norm() == Approx(0.0).margin(1e-9));
    condition = (imm.initializator.filterBal->getPredictInfo().Pe - expectedPredictStateCovarianceBalreentry).norm() == Approx(0.0).margin(1e-9);
    if(!condition){
        PRINTM(imm.initializator.filterBal->getPredictInfo().Pe);
        PRINTM(expectedPredictStateCovarianceBalreentry);
    }


    Eigen::MatrixXd expectedXpred (9,1);
    expectedXpred << -208412.351684929,475.553922395592,-0.0607498364397687,155930.573414262,-469.708173062171,0.0457523856592994,26844.908135732,2130.98411236917,0.00795357406674756;
    Eigen::MatrixXd expectedPpred (9,9);
    expectedPpred << 423729.253149812,839997.112987463,-38.5291965180739,534862.07566649,1070100.66550604,30.160003679966,233601.044174825,500619.051557954,5.15347121453614,
                    839997.112987463,2186004.76785978,-77.0844490124841,1064171.29300814,2788202.87555395,66.253573420642,502962.932904809,1190589.220522,11.3360725224568,
                    -38.5291965180739,-77.0844490124841,100.877011633307,-61.0016405284388,-137.931685608594,-0.0284633901493936,31.9189590268445,71.171067258101,-0.00447628515134929,
                    534862.07566649,1064171.29300814,-61.0016405284388,725424.765789859,1444036.08818491,47.4148178756333,61268.9607798526,152926.159331743,7.66195144473268,
                    1070100.66550604,2788202.87555395,-137.931685608594,1444036.08818491,3783781.46473397,114.691109280077,150649.01599803,287055.465453587,17.3156826773172,
                    30.160003679966,66.253573420642,-0.0284633901493936,47.4148178756333,114.691109280077,100.861578707748,-23.5369458708332,-52.4784479333816,0.00342830734514459,
                    233601.044174825,502962.932904809,31.9189590268445,61268.9607798526,150649.01599803,-23.5369458708332,1475960.04382186,3021780.06317346,-1.60597027421826,
                    500619.051557954,1190589.220522,71.171067258101,152926.159331743,287055.465453587,-52.4784479333816,3021780.06317346,7603938.148565,3.94222461491343,
                    5.15347121453614,11.3360725224568,-0.00447628515134929,7.66195144473268,17.3156826773172,0.00342830734514459,-1.60597027421826,3.94222461491343,100.840352169794;

    // CHECK((imm.getPredictInfo().Xe - expectedXpred).norm() == Approx(0.0).margin(1e-9));
    // condition = (imm.getPredictInfo().Xe - expectedXpred).norm() == Approx(0.0).margin(1e-9);
    // if(!condition){
    //     PRINTM(imm.getPredictInfo().Xe);
    //     PRINTM(expectedXpred);
    // }

    // CHECK((imm.getPredictInfo().Pe - expectedPpred).norm() == Approx(0.0).margin(1e-9));
    // condition = (imm.getPredictInfo().Pe - expectedPpred).norm() == Approx(0.0).margin(1e-9);
    // if(!condition){
    //     PRINTM(imm.getPredictInfo().Pe);
    //     PRINTM(expectedPpred);
    // }


    BENCHMARK("Imm_predict"){
    // auto pred = imm.predict(dt);
    };

}
    

   