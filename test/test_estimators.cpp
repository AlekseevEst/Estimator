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
    X << -1255.5, 0.0, -18.039, 0.0, 18.0, 0.0; // тест нерабочий!!!
    Z << -1258.5, 3.119, 1.885;                 // тест нерабочий!!!
    Q << 0.5,0.0,0.0,
        0.0,0.5,0.0,
        0.0,0.0,0.5;
    R << 1.0,0.0,0.0,
        0.0,1e-4,0.0,
        0.0,0.0,1e-4;

    double t = 6.0;
    double k = 1.0;


    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSph> ukf (X,Z,t,Q,R,k);
    CHECK(ukf.predict() == X);
    CHECK(ukf.correct(Z) == X);

    BENCHMARK("ESTIMATOR"){
        ukf.predict();
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
