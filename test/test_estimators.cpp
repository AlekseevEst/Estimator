#include <catch2/catch.hpp>
#include "ukf.h"
using namespace Catch::Benchmark;

TEST_CASE("filtering_ukf_CV")
{

    Eigen::MatrixXd X(6, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(3, 3);
    Eigen::MatrixXd R(3, 3);
    X << 1200.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Z << 2400.0, 1e-02, 1e-02;                 

    Q << 0.5,0.0,0.0,
        0.0,0.5,0.0,
        0.0,0.0,0.5;
        
    R << 1.0,0.0,0.0,
        0.0,1e-4,0.0,
        0.0,0.0,1e-4;


    double t = 6.0;
    double k = 1.0;

    Eigen::MatrixXd expectedCorrectState(6,1);
    Eigen::MatrixXd expectedPredState(6,1);
    
    expectedCorrectState << 2389.32, 198.165, 11.6005, 1.93207, 11.6005, 1.93207;
    expectedPredState << 1200.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph> ukf (X,t,Q,R,k);
    CHECK(ukf.predict().isApprox(expectedPredState,0.0001));
    CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.0001));

    BENCHMARK("STEP"){
        ukf.predict();
        ukf.correct(Z);
    };
}

TEST_CASE("filtering_ukf_CT")
{

    Eigen::MatrixXd X(7, 1);
    Eigen::MatrixXd Z(4, 1);
    Eigen::MatrixXd Q(4, 4);
    Eigen::MatrixXd R(4, 4);
    X << 1200.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Z << 2400.0, 1e-02, 1e-02, 200.0;                 

    Q << 0.5,0.0,0.0,0.0,
        0.0,0.5,0.0,0.0,
        0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0001;
        
    R << 1.0, 0.0, 0.0, 0.0,
        0.0, 1e-4, 0.0, 0.0,
        0.0, 0.0, 1e-4, 0.0,
        0.0, 0.0, 0.0, 10.0;


    double t = 6.0;
    double k = 1.0;

    Eigen::MatrixXd expectedCorrectState(7,1);
    Eigen::MatrixXd expectedPredState(7,1);
    
    expectedCorrectState << 2400.0, 198.165, 11.6005, 1.93207, 11.6005, 1.93207,0.0;
    expectedPredState << 1200.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSph> ukf (X,t,Q,R,k);
    CHECK(ukf.predict().isApprox(expectedPredState,0.0001));
    CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.0001));


    BENCHMARK("STEP"){
        ukf.predict();
        ukf.correct(Z);
    };
}
