#pragma once

#include <array>
#include <cmath>
#include <vector>
#include <iostream>
#include "utils.h"
#include "models.h"
#include "ukfMath.h"

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>

struct UnscentedKalmanfilter
{
private:
    UnscentedKalmanFilterMath<M> UKfilterMath;
    double T;
    std::vector<double> weightVectors;
    M extrapolatedStateSigmaVectors;
    StateFunc<M> stateFunc;
    MeasurementFunc<M> measFunc;

public:
    Predict<M> predictStruct;
    Correct<M> correctStruct;

    M predict();
    M correct(const M &Z);

    UnscentedKalmanfilter(const M& X, double t, const M& procNoise, const M& measNoiseMatRadian ,double koef): 
                                                                                                            UKfilterMath(measNoiseMatRadian, procNoise, koef)    
                                                                                                            {   
                                                                                                                T = t;
                                                                                                                correctStruct.X = X;
                                                                                                                correctStruct.P.resize(X.rows(),X.rows());
                                                                                                                correctStruct.P.setZero();
                                                                                                            }
};

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>

M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc>::predict()
{
    correctStruct.P = UKfilterMath.make_P_cart(correctStruct.P, correctStruct.X);
    M sqrtOfP = UKfilterMath.sqrt_matrix_p(correctStruct.P);
    M sigmaVectors = UKfilterMath.doSigmaVectors(correctStruct.X,sqrtOfP);
    weightVectors = UKfilterMath.calculationVectorWeights();
    extrapolatedStateSigmaVectors = stateFunc(sigmaVectors, T);
    predictStruct.Xe = UKfilterMath.doExtrapolatedStateVector(extrapolatedStateSigmaVectors, weightVectors);
    predictStruct.Pe = UKfilterMath.doCovMatExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictStruct.Xe, weightVectors, T);

    return predictStruct.Xe;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>
M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc>::correct(const M& Z)
{   
    M extrapolatedMeasSigmaVectors = measFunc(extrapolatedStateSigmaVectors, Z);
    predictStruct.Ze = UKfilterMath.doExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, weightVectors);
    predictStruct.Se = UKfilterMath.doCovMatExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictStruct.Ze, weightVectors);
    predictStruct.K = UKfilterMath.calcGainFilter(extrapolatedStateSigmaVectors, predictStruct.Xe, extrapolatedMeasSigmaVectors,predictStruct.Ze, predictStruct.Se, weightVectors);
    correctStruct.X = UKfilterMath.correctState(predictStruct.Xe,Z, predictStruct.Ze, predictStruct.K);
    correctStruct.P = UKfilterMath.correctCov(predictStruct.Pe, predictStruct.K, predictStruct.Se);

    return correctStruct.X;
}
