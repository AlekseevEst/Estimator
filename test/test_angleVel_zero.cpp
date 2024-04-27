#include <catch2/catch.hpp>
#include "ukf.h"
using namespace Catch::Benchmark;

TEST_CASE("angleVel_zero")
{
   Eigen::MatrixXd X(7, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(4, 4);
    Eigen::MatrixXd R(3, 3);
    X << 50.0, 200.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    Z << 100.0, 1e-02, 1e-02;                 

    Q << 0.5,0.0,0.0,0.0,
        0.0,0.5,0.0,0.0,
        0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0;
        
    R << 1.0, 0.0, 0.0,
        0.0, 1e-4, 0.0,
        0.0, 0.0, 1e-4;


    double t = 0.25;
    double k = 1.0;

    Eigen::MatrixXd expectedPredState(7,1);
    Eigen::MatrixXd expectedCorrectState(7,1);

    expectedPredState << 100.0, 200.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    expectedCorrectState << 99.945, 199.811, 0.863, 3.450, 0.863, 3.450, 0.0;
    Points p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = 0;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSph> ukf (X,t,Q,R,p);

    CHECK(ukf.predict().isApprox(expectedPredState,0.0001));
    CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.0001));

    BENCHMARK("angleVel_zero")
    {
        ukf.predict();
        ukf.correct(Z);
    };
}
