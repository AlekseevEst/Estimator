#include <catch2/catch.hpp>
#include "ukf.h"
using namespace Catch::Benchmark;


TEST_CASE("filtering_ukf_CT")
{

    Eigen::MatrixXd X(7, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(4, 4);
    Eigen::MatrixXd R(3, 3);
    X << 9922.94637347, 0.0, 19978.67752417, 0.0, 9955.79978476, 0.0, 0.0;
    Z << 24544.6547, 63.6329241, 24.192866;                 

    Q << 10.0,0.0,0.0,0.0,
        0.0,10.0,0.0,0.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1e-7;
        
    R << 10000.0, 0.0, 0.0,
        0.0, pow((0.1/3),2), 0.0,
        0.0, 0.0, pow((0.1/3),2);


    double t = 0.25;

    Eigen::MatrixXd expectedCorrectState(7,1);
    Eigen::MatrixXd expectedPredState(7,1);
    
    expectedPredState << 9922.95, 0.0, 19978.7, 0.0, 9955.8, 0.0, 0.0;
    expectedCorrectState << 9933.22, 0.0, 20019.3, 0.0, 10035.2, 224.576, 0.0;
    
    Points p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = -4;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSphCVCT> ukf (X,Q,R,p);
    // std::cout<<ukf.predict();
    // std::cout<<ukf.correct(Z);
    CHECK(ukf.predict(t).isApprox(expectedPredState,0.00001));
    CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.00001));

    BENCHMARK("STEP"){
        ukf.predict(t);
        ukf.correct(Z);
    };
}
// TEST_CASE("filtering_ukf_CV")
// {

//     Eigen::MatrixXd X(6, 1);
//     Eigen::MatrixXd Z(3, 1);
//     Eigen::MatrixXd Q(3, 3);
//     Eigen::MatrixXd R(3, 3);
//     X << 1200.065, 0.0, 105.6025, 0.0, 50.3206, 0.0;
//     Z << 2401.302, 1e-02, 1e-02;                 

//     Q << 0.5,0.0,0.0,
//         0.0,0.5,0.0,
//         0.0,0.0,0.5;
        
//     R << 1.0,0.0,0.0,
//         0.0,1e-4,0.0,
//         0.0,0.0,1e-4;


//     double t = 6.0;
//     double k = 1.0;

//     Eigen::MatrixXd expectedCorrectState(6,1);
//     Eigen::MatrixXd expectedPredState(6,1);
    
//     expectedCorrectState << 2390.23, 198.305, 119.923, 2.39244, 63.2484, 2.15664;
//     expectedPredState << 1200.07, 0.0, 105.602, 0.0, 50.3206, 0.0;

//     Points p;
//     p.alpha = 1e-3;
//     p.beta = 2;
//     p.kappa = 0;

//     UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph> ukf (X,t,Q,R,p);
//     CHECK(ukf.predict().isApprox(expectedPredState,0.0001));
//     CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.0001));

//     BENCHMARK("STEP"){
//         ukf.predict();
//         ukf.correct(Z);
//     };
// }
