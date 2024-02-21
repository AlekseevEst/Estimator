#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "ukf.h"
#include "ekf.h"
using namespace Catch::Benchmark;

TEST_CASE("filtering_ukf")
{
    // UnscentKalmanfilter<> uf;
    // Eigen::MatrixXd X(6, 1);
    // Eigen::MatrixXd Z(3, 1);
    // X << -1255.5, 0.0, -18.039, 0.0, 18.0, 0.0; // тест нерабочий!!!
    // Z << -1258.5, 3.119, 1.885;                 // тест нерабочий!!!
    // CHECK(uf.predictUkf(X) == X);
    // CHECK(uf.correctUkf(Z) == X);

    // CHECK();
    BENCHMARK("ESTIMATOR"){
        // uf.predictUkf(X);
    };
}

TEST_CASE("filtering_ekf")
{
    // ExtendedKalmanFilter<> ekf;
    // Eigen::MatrixXd X(6, 1);
    // Eigen::MatrixXd Z(3, 1);
    // X << -3000.5, 0.0, -50.039, 0.0, 80.0, 0.0; // тест нерабочий!!!
    // Z << -3020.5, 3.0, 1.7;                 // тест нерабочий!!!
    // CHECK(ekf.predictEkf(X) == X);
    // CHECK(ekf.correctEkf(Z) == X);

    // CHECK();
    BENCHMARK("ESTIMATOR"){
        // ekf.predictEkf(X);
    };
}
