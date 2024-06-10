#include <catch2/catch.hpp>
#include "track.h"
using namespace Catch::Benchmark;

TEST_CASE("Track")
{

    Eigen::MatrixXd Z0(3, 1);
    Eigen::MatrixXd Z1(3, 1);

    Z0 << 24544.6547, 63.6329241, 24.192866;                 
    Z1 << 24844.6547, 65.6329241, 28.192866;    

    Detection<Eigen::MatrixXd> d0;
    Detection<Eigen::MatrixXd> d1;
    d0.timePoint = 0.25;
    d0.point = Z0;
    d1.timePoint = 0.5;
    d1.point = Z1;
    Eigen::MatrixXd expectedCorrectState(9,1);
    Eigen::MatrixXd expectedPredState(9,1);
    
    // expectedPredState << 9922.95, 0.0, 19978.7, 0.0, 9955.8, 0.0, 0.0;
    // expectedCorrectState << 9933.22, 0.0, 20019.3, 0.0, 10035.2, 224.576, 0.0;
    

    Track <Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>, 
                                                InitUKFStateModelCAMeasureModelSph<Eigen::MatrixXd,FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>> track (d0);
  
    Eigen::MatrixXd correct = track.step(d1);
    PRINTM(correct);
    // std::cout<<ukf.predict();
    // std::cout<<ukf.correct(Z);
    // CHECK(ukf.predict(t).isApprox(expectedPredState,0.00001));
    // CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.00001));

    BENCHMARK("STEP"){

    };
}