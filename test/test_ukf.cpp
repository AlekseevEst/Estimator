// #include <iostream>
// #include <catch2/catch.hpp>
// #include "imm.h"
// #include "ukf.h"
// #include "models.h"
// // #include "ifilter.h"
// // #include "initFilters.h"
// using namespace Catch::Benchmark;


// struct InitUnscentedKalmanFilterCVtest
// {
//     template <class M,
//               class... TypeArgs>
//     void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCVtest, TypeArgs...> &filter,
//                     const M &meas,
//                     const M &measNoise)
//     {
      
//         filter.correctInfo.X << 20. ,2. ,20. ,2., 20.,2.;



//         filter.correctInfo.P << 1., 0, 0, 0, 0, 0,
//                                 0, 1. ,0 ,0 ,0 ,0,
//                                 0, 0, 1. ,0 , 0 ,0,
//                                 0, 0, 0 ,1. , 0 ,0,
//                                 0, 0, 0 ,0 , 1. ,0,
//                                 0, 0, 0 ,0 , 0 ,1.;

//             filter.Q.resize(3,3);
//             filter.Q << 1.0, 0.0, 0.0,
//                         0.0, 1.0, 0.0,
//                         0.0, 0.0, 1.0;

//         filter.R = measNoise;

//         filter.paramsSigmaPoints.alpha = 1e-3;
//         filter.paramsSigmaPoints.beta = 2.0;
//         filter.paramsSigmaPoints.kappa = 0.0;
//     }
// };



// TEST_CASE("imm_Predict")
// {
//     Eigen::MatrixXd meas(3,1);
//     Eigen::MatrixXd measNoise(3,3);
//     measNoise << 1.,0,0,
//                  0,1.,0,
//                  0,0,1.;
 


//     UnscentedKalmanFilter<Eigen::MatrixXd, InitUnscentedKalmanFilterCVtest, FuncConstVel<Eigen::MatrixXd>, FuncMeasSph<Eigen::MatrixXd>, FuncControlMatrix_XvXYvYZvZ<Eigen::MatrixXd>> ukf;
//     ukf.Initialization(meas, measNoise);
//     double dt = 0.2;



//     Eigen::MatrixXd expectedXpred(6,1);
//     expectedXpred << 20.4, 2., 20.4, 2., 20.4,  2.;

//     Eigen::MatrixXd expectedPpred(6,6);
//     expectedPpred <<1.04040000e+00,  2.04000000e-01,  2.43642889e-17,  4.56864242e-20, 4.96478152e-17,  5.21723947e-18,
//                     2.04000000e-01,  1.04000000e+00, -1.10352854e-18, -2.15397832e-21, -2.24862598e-18, -2.36398293e-19,
//                     2.43642201e-17, -1.10351257e-18,  1.04040000e+00,  2.04000000e-01, 4.96478152e-17,  5.21723947e-18,
//                     4.56960271e-20, -2.15381345e-21,  2.04000000e-01,  1.04000000e+00, 9.30996657e-20,  9.68722481e-21,
//                     4.96475978e-17, -2.24854928e-18,  4.96478095e-17,  9.31268867e-20, 1.04040000e+00,  2.04000000e-01,
//                     5.21724996e-18, -2.36401479e-19,  5.21724996e-18,  9.69290242e-21, 2.04000000e-01,  1.04000000e+00;


//     Eigen::MatrixXd expectedXcor(6,1);
//     expectedXcor << 18.60054543,  1.64716577, 18.60054543,  1.64716577, 23.99389129,  2.70468457;
//     Eigen::MatrixXd expectedPcor(6,6);
//     expectedPcor << 0.3185912,  0.06246886, 0.11473142, 0.02249636, 0.07718949, 0.01513519,
//                     0.06246886, 1.0122488,  0.02249636, 0.00441105, 0.01513519, 0.00296769,
//                     0.11473142, 0.02249636, 0.3185912,  0.06246886, 0.07718949, 0.01513519,
//                     0.02249636, 0.00441105, 0.06246886, 1.0122488,  0.01513519, 0.00296769,
//                     0.07718949, 0.01513519, 0.07718949, 0.01513519, 0.35564801, 0.0697349, 
//                     0.01513519, 0.00296769, 0.01513519, 0.00296769, 0.0697349, 1.01367351;                    


//     auto pred = ukf.predict(dt);

//     CHECK((pred.first - expectedXpred).norm() == Approx(0.0).margin(1e-2));
//     bool condition = (pred.first - expectedXpred).norm() == Approx(0.0).margin(1e-2);
//     if(!condition){
//         PRINTM(pred.first.format(customFormat));
//         PRINTM(expectedXpred.format(customFormat));
//     }

//     CHECK((pred.second - expectedPpred).norm() == Approx(0.0).margin(1e-2));
//     condition = (pred.second - expectedPpred).norm() == Approx(0.0).margin(1e-2);
//     if(!condition){
//         PRINTM(pred.second.format(customFormat));
//         PRINTM(expectedPpred.format(customFormat));
//     }

//     Eigen::MatrixXd Z (3,1);
//     Z << 35.355, 45., 45.; 
    
//     auto cor = ukf.correct(Z);

//     CHECK((cor.first - expectedXcor).norm() == Approx(0.0).margin(1e-1));
//     condition = (cor.first - expectedXcor).norm() == Approx(0.0).margin(1e-1);
//     if(!condition){
//         PRINTM(cor.first.format(customFormat));
//         PRINTM(expectedXcor.format(customFormat));
//     }

//     CHECK((cor.second - expectedPcor).norm() == Approx(0.0).margin(1e-1));
//     condition = (cor.second - expectedPcor).norm() == Approx(0.0).margin(1e-1);
//     if(!condition){
//         PRINTM(cor.second.format(customFormat));
//         PRINTM(expectedPcor.format(customFormat));
//     }


//     double expectedlikelihood = 3.0182580264650036e-08;
//     double likelihood = ukf.likelihood(Z);
//     CHECK(likelihood - expectedlikelihood == Approx(0.0).margin(1e-10));
//     condition = likelihood - expectedlikelihood == Approx(0.0).margin(1e-10);
//     if(!condition){
//         std::cout<< std::fixed<<std::setprecision(20)<< likelihood<<std::endl;;
//         std::cout<< std::fixed<<std::setprecision(20)<< expectedlikelihood<<std::endl;;
//     }

//     BENCHMARK("Imm_predict"){
//     auto pred = ukf.predict(dt);
//     auto cor = ukf.correct(Z);
//     };

// }
    

   