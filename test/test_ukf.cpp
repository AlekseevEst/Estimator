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

    Points p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = 0;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph> ukf (X,t,Q,R,p);
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
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(4, 4);
    Eigen::MatrixXd R(3, 3);
    X << 49.9, 200.0, 0.0, 0.0, 0.0, 0.0, 0.098;
    Z << 99.9, 1e-01, 1e-01;                 

    Q << 0.5,0.0,0.0,0.0,
        0.0,0.5,0.0,0.0,
        0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0;
        
    R << 1.0, 0.0, 0.0,
        0.0, 1e-4, 0.0,
        0.0, 0.0, 1e-4;


    double t = 0.25;
    double k = 1.0;

    Eigen::MatrixXd expectedCorrectState(7,1);
    Eigen::MatrixXd expectedPredState(7,1);
    
    expectedPredState << 99.895, 199.94, 0.612, 4.895, 0.0, 0.0, 0.098;
    expectedCorrectState << 99.792, 199.189, 8.708, 37.255, 8.625, 34.477, 0.098;
    
    Points p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = 0;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSph> ukf (X,t,Q,R,p);
    CHECK(ukf.predict().isApprox(expectedPredState,0.0001));
    CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.0001));


    BENCHMARK("STEP"){
        ukf.predict();
        ukf.correct(Z);
    };
}
