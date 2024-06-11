#include <catch2/catch.hpp>
#include "ukf.h"
using namespace Catch::Benchmark;


TEST_CASE("filtering_ukf_CT")
{

    Eigen::MatrixXd X(7, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(4, 4);
    Eigen::MatrixXd R(3, 3);
    X << 9922.94637347, 0.0, 19978.67752417, 0.0, 9955.79978476, 0.0, 0.0;
    Z << 24544.6547, 63.6329241, 24.192866;                 

    Q << 10.0,0.0,0.0,0.0,
        0.0,10.0,0.0,0.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1e-7;
        
    R << 10000.0, 0.0, 0.0,
        0.0, pow((0.1/3),2), 0.0,
        0.0, 0.0, pow((0.1/3),2);


    double t = 0.25;

    Eigen::MatrixXd expectedCorrectState(7,1);
    Eigen::MatrixXd expectedPredState(7,1);
    
    expectedPredState << 9922.95, 0.0, 19978.7, 0.0, 9955.8, 0.0, 0.0;
    expectedCorrectState << 9933.22, 0.0, 20019.3, 0.0, 10035.2, 224.576, 0.0;
    
    ParamSigmaPoints p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = -4;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW> ukf (X,Q,R,p);
    // std::cout<<ukf.predict();
    // std::cout<<ukf.correct(Z);
    CHECK(ukf.predict(t).isApprox(expectedPredState,0.00001));
    CHECK(ukf.correct(Z).isApprox(expectedCorrectState,0.00001));

    BENCHMARK("STEP"){
        ukf.predict(t);
        ukf.correct(Z);
    };
}
TEST_CASE("filtering_ukf_CV_with_Vr")
{

    Eigen::MatrixXd X(6, 1);
    Eigen::MatrixXd Z(4, 1);
    Eigen::MatrixXd Q(3, 3);
    Eigen::MatrixXd R(4, 4);
    X << 130000.6547, 0.0, 0.6329241, 0.0, 0.192866, 0.0;
    Z << 131200.1248, 0.5256113, 0.182648, 189.0;                 

    Q << 0.00001,        0.0,        0.0,
             0.0,    0.00001,        0.0,
             0.0,        0.0,        0.00001;
 
        
    R << 10000.0,         0.0,                0.0,          0.0,
             0.0,    pow((0.1/3),2),          0.0,          0.0,
             0.0,         0.0,           pow((0.1/3),2),    0.0,
             0.0,         0.0,                0.0,       pow(5.0 ,2);

    double t = 6.0;

    Eigen::MatrixXd expectedCorrectState(6,1);
    Eigen::MatrixXd expectedPredState(6,1);
    
    expectedPredState << 130001, 0.0, 0.6329, 0.0, 0.192866, 0.0;
    expectedCorrectState << 130324, 2.31721, 1192.58, 198.657, 414.418, 69.0373;
    
    ParamSigmaPoints p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = -3;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ> ukf (X,Q,R,p);
    // std::cout<<ukf.predict(t);
    // std::cout<<ukf.correct(Z);
    Eigen::MatrixXd predict = ukf.predict(t);
    Eigen::MatrixXd correct = ukf.correct(Z);

    CHECK(predict.isApprox(expectedPredState,0.00001));
    CHECK(correct.isApprox(expectedCorrectState,0.00001));

    BENCHMARK("STEP"){
        ukf.predict(t);
        ukf.correct(Z);
    };
}

TEST_CASE("filtering_ukf_CA")
{

    Eigen::MatrixXd X(9, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(3, 3);
    Eigen::MatrixXd R(3, 3);
    X <<  22506.5388, 0.0, 0.0,  8.01622656e-03, 0.0, 0.0,  6.29699009e+01, 0.0, 0.0;
    Z << 22408.2757, 1.55981886e-02, 6.24263741e+01;                 

    Q << 10.0,  0.0,    0.0,
          0.0,  10.0,   0.0,
         0.0,   0.0,    10.0;

        
    R << 10000.0,           0.0,            0.0,
           0.0,        pow((0.1/3),2),      0.0,
           0.0,             0.0,        pow((0.1/3),2);


    double t = 1.0;

    Eigen::MatrixXd expectedCorrectState(9,1);
    Eigen::MatrixXd expectedPredState(9,1);
    
    expectedPredState << 22506.5, 0.0, 0.0, 0.00801623, 0.0, 0.0,  62.9699, 0.0, 0.0;
    expectedCorrectState << 22352.1, -125.877, -3.81446, 6.10151, 6.16124, 0.186704, 24419.3, 24627.1, 746.276;
    
    ParamSigmaPoints p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = -6;

    UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ> ukf (X,Q,R,p);
    // std::cout<<ukf.predict(t);
    // std::cout<<ukf.correct(Z);
    Eigen::MatrixXd predict = ukf.predict(t);
    Eigen::MatrixXd correct = ukf.correct(Z);

    CHECK(predict.isApprox(expectedPredState,0.00001));
    CHECK(correct.isApprox(expectedCorrectState,0.00001));

    BENCHMARK("STEP"){
        ukf.predict(t);
        ukf.correct(Z);
    };
}