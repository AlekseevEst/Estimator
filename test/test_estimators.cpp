#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "ukf.h"
#include "ekf.h"
using namespace Catch::Benchmark;

TEST_CASE("filtering_ukf")
{


    Eigen::MatrixXd X(6, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(3, 3);
    Eigen::MatrixXd R(3, 3);
    X << 1700.24, 0.0, 17.01, 0.0, 17.02, 0.0;
    Z << 2888.5, 1e-02, 1e-02;                 

    Q << 0.5,0.0,0.0,
        0.0,0.5,0.0,
        0.0,0.0,0.5;
    R << 1.0,0.0,0.0,
        0.0,1e-4,0.0,
        0.0,0.0,1e-4;

    double t = 6.0;
    double k = 1.0;

    Eigen::MatrixXd x_check(6,1);
    x_check << 2.88476682e+03, 1.97394641e+02, 2.90383652e+01, 2.00563475e+00, 2.90209862e+01, 2.00259655e+00;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph> ukf (X,Z,t,Q,R,k);
    ukf.predict();
    CHECK(ukf.correct(Z).isApprox(x_check, 1.0));

    BENCHMARK("ESTIMATOR"){
        ukf.correct(Z);
    };
}

// TEST_CASE("filtering_ekf")
// {
    // ExtendedKalmanFilter<> ekf;
    // Eigen::MatrixXd X(6, 1);
    // Eigen::MatrixXd Z(3, 1);
    // X << -3000.5, 0.0, -50.039, 0.0, 80.0, 0.0; // тест нерабочий!!!
    // Z << -3020.5, 3.0, 1.7;                 // тест нерабочий!!!
    // CHECK(ekf.predictEkf(X) == X);
    // CHECK(ekf.correctEkf(Z) == X);

    // CHECK();
    // BENCHMARK("ESTIMATOR"){
        // ekf.predictEkf(X);
    // };
// }
