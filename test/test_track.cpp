#include <catch2/catch.hpp>
#include "track.h"
using namespace Catch::Benchmark;

TEST_CASE("Track_CA_MeasSph")
{

    Eigen::MatrixXd Z0(3, 1);
    Eigen::MatrixXd Z1(3, 1);

    Z0 << 22506.5388, 8.01622656e-03, 6.29699009e+01;                 
    Z1 << 22408.2757, 1.55981886e-02, 6.24263741e+01;    

    Detection<Eigen::MatrixXd> d0;
    Detection<Eigen::MatrixXd> d1;
    d0.timePoint = 1.0;
    d0.point = Z0;
    d1.timePoint = 2.0;
    d1.point = Z1;
    Eigen::MatrixXd expectedCorrectState(9,1);
    Eigen::MatrixXd expectedPredState(9,1);
    
    expectedCorrectState << 10381, 124.44, 3.77091, 2.8046, 1.39348, 0.04222, 19882.9, -167.619, -5.07937;
    

    Track <Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>, 
                                                InitUKFStateModelCAMeasureModelSph<Eigen::MatrixXd,FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>> track (d0);
  
    Eigen::MatrixXd correct = track.step(d1);

    CHECK(correct.isApprox(expectedCorrectState,0.00001));
    BENCHMARK("STEP"){

    };
}

TEST_CASE("Track_CT_MeasSph")
{

    Eigen::MatrixXd Z0(3, 1);
    Eigen::MatrixXd Z1(3, 1);

    Z0 << 24613.0220, 63.5438, 24.0914;                 
    Z1 << 24352.4686, 63.6769, 24.1093;    

    Detection<Eigen::MatrixXd> d0;
    Detection<Eigen::MatrixXd> d1;
    d0.timePoint = 0.25;
    d0.point = Z0;
    d1.timePoint = 0.5;
    d1.point = Z1;
    Eigen::MatrixXd expectedCorrectState(7,1);

    expectedCorrectState << 9933.23, 0.0, 20019.8, 0.0, 9999.65, 19.4815, 0.0;
    

    Track <Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>, 
                                                InitUKFStateModelCTMeasureModelSph<Eigen::MatrixXd,FuncConstTurn, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>> track (d0);
  
    Eigen::MatrixXd correct = track.step(d1);

    CHECK(correct.isApprox(expectedCorrectState,0.00001));
    BENCHMARK("STEP"){

    };
}

TEST_CASE("Track_CV_MeasSph")
{

    Eigen::MatrixXd Z0(3, 1);
    Eigen::MatrixXd Z1(3, 1);

    Z0 << 130000.6547, 0.6329241, 0.192866;                 
    Z1 << 131200.1248, 0.5256113, 0.182648;      

    Detection<Eigen::MatrixXd> d0;
    Detection<Eigen::MatrixXd> d1;
    d0.timePoint = 6.0;
    d0.point = Z0;
    d1.timePoint = 12.0;
    d1.point = Z1;
    Eigen::MatrixXd expectedCorrectState(6,1);
    
    expectedCorrectState << 131175.0, 195.757, 1206.55, -38.1009, 418.499, -3.17293;
    

    Track <Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>, 
                                                InitUKFStateModelCVMeasureModelSph<Eigen::MatrixXd,FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>> track (d0);
  
    Eigen::MatrixXd correct = track.step(d1);
    CHECK(correct.isApprox(expectedCorrectState,0.00001));

    BENCHMARK("STEP"){

    };
}

TEST_CASE("Track_CV_with_Vr_MeasSph")
{

    Eigen::MatrixXd Z0(4, 1);
    Eigen::MatrixXd Z1(4, 1);

    Z0 << 130000.6547, 0.6329241, 0.192866, 190.0;                 
    Z1 << 131200.1248, 0.5256113, 0.182648, 189.0;    

    Detection<Eigen::MatrixXd> d0;
    Detection<Eigen::MatrixXd> d1;
    d0.timePoint = 6.0;
    d0.point = Z0;
    d1.timePoint = 12.0;
    d1.point = Z1;
    Eigen::MatrixXd expectedCorrectState(6,1);

    expectedCorrectState << 130318, 2.778, 1196.13, -40.5525, 416.458, -3.6972;
    

    Track <Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>, 
                                                InitUKFStateModelCVMeasureModelSph<Eigen::MatrixXd,FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>> track (d0);
  
    Eigen::MatrixXd correct = track.step(d1);

    CHECK(correct.isApprox(expectedCorrectState,0.00001));

    BENCHMARK("STEP"){

    };
}