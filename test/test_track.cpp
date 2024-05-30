#include <catch2/catch.hpp>
#include "track.h"
using namespace Catch::Benchmark;

TEST_CASE("Track")
{

    // Eigen::MatrixXd X(7, 1);
    // Eigen::MatrixXd Z(3, 1);
    // Eigen::MatrixXd Q(4, 4);
    // Eigen::MatrixXd R(3, 3);
    // X << 9922.94637347, 0.0, 19978.67752417, 0.0, 9955.79978476, 0.0, 0.0;
    // Z << 24544.6547, 63.6329241, 24.192866;                 

    // Q << 10.0,0.0,0.0,0.0,
    //     0.0,10.0,0.0,0.0,
    //     0.0,0.0,1.0,0.0,
    //     0.0,0.0,0.0,1e-7;
        
    // R << 10000.0, 0.0, 0.0,
    //     0.0, pow((0.1/3),2), 0.0,
    //     0.0, 0.0, pow((0.1/3),2);


    // double t = 0.25;

    // Eigen::MatrixXd expectedCorrectState(7,1);
    // Eigen::MatrixXd expectedPredState(7,1);
    
    // expectedPredState << 9922.95, 0.0, 19978.7, 0.0, 9955.8, 0.0, 0.0;
    // expectedCorrectState << 9933.22, 0.0, 20019.3, 0.0, 10035.2, 224.576, 0.0;
    
    // Points p;
    // p.alpha = 1e-3;
    // p.beta = 2;
    // p.kappa = -4;

    // UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSphCVCT> ukf (X,Q,R,p);
    // // std::cout<<ukf.predict();
    // // std::cout<<ukf.correct(Z);
    // CHECK(ukf.predict(t).isApprox(expectedPredState,0.00001));
    // CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.00001));

    // BENCHMARK("STEP"){

    // };
}