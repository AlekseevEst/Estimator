#include <catch2/catch.hpp>
#include "ukf.h"
using namespace Catch::Benchmark;


TEST_CASE("angle")
{

    double diff = Utils<Eigen::MatrixXd>::ComputeAngleDifference(-1.0 * (M_PI/180.0), 181.0 * (M_PI/180.0)) * (180.0/M_PI);
    std::cout<<diff;
    double a = 178.000000;

    CHECK(a == Approx(diff).epsilon(0.00001));

    BENCHMARK("angle"){

    };
}

