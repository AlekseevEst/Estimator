#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iomanip>
#include "ukf.h"
#include "Eigen/Dense"
using namespace Catch::Benchmark;

TEST_CASE("Predict")
{   unscent_filter uf;
    MatrixXd X(6,1);
    MatrixXd Z(3,1);
    X <<-1255.5,0.0, -18.039, 0.0, 18.0, 0.0; // тест нерабочий!!!
    Z << -1258.5,3.119, 1.885;                // тест нерабочий!!!
    CHECK(uf.predictUkf(X) == X);
    CHECK(uf.correctUkf(Z) == X);

    // CHECK();
    BENCHMARK("ESTIMATOR")
    {
        // uf.predictUkf(X);
    };
}

