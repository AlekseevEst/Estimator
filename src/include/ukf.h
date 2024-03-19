#pragma once

#include <array>
#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <chrono>
#include <boost/math/distributions/chi_squared.hpp>
#include "utils.h"
#include "models.h"
#include "ukfMath.h"

template <class M,
                template <typename> class StateFunc,
                template <typename> class MeasurementFunc>

struct UnscentedKalmanfilter
{
    private:

    double T;

    UnscentedKalmanFilterMath<M>  UKfilterMath;
    std::vector<double> weightVectors;
    std::vector<M> ExtrapolatedSigmaVectors;
    StateFunc<M> stateFunc;
    MeasurementFunc<M> measFunc;
    
    public:

    M predict();
    M correct(const M& Z);

    UnscentedKalmanfilter(const M& X, double t, const M& procNoise, const M& measNoiseMatRadian ,double koef): 
                                                                                                            T(t),UKfilterMath(X,measNoiseMatRadian, procNoise, t, koef){}
};

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>

M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc>::predict()
{
    UKfilterMath.make_P_cart();
    try
    {
        UKfilterMath.sqrt_matrix_p();
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << e.what() << '\n';
    }
    std::vector<M> sigmaVectors = UKfilterMath.doSigmaVectors();
    weightVectors = UKfilterMath.calculationVectorWeights();
    ExtrapolatedSigmaVectors = stateFunc(sigmaVectors, T);
    UKfilterMath.doExtrapolatedStateVector(ExtrapolatedSigmaVectors, weightVectors);
    UKfilterMath.doCovMatExtrapolatedStateVector(ExtrapolatedSigmaVectors, weightVectors);
    
    return UKfilterMath.predictStruct.Xe;
}

template <class M,
          template <typename> class StateFunc,
          template <typename> class MeasurementFunc>
M UnscentedKalmanfilter<M, StateFunc, MeasurementFunc>::correct(const M& Z)
{   
    std::vector<M> ExtrapolatedSigmaVectorsSph = measFunc(ExtrapolatedSigmaVectors, Z);
    UKfilterMath.doExtrapolatedStateVectorSph(ExtrapolatedSigmaVectorsSph, weightVectors);
    UKfilterMath.doCovMatExtrapolatedStateVectorSph(ExtrapolatedSigmaVectorsSph, weightVectors);
    UKfilterMath.calcGainFilter(ExtrapolatedSigmaVectors, ExtrapolatedSigmaVectorsSph, weightVectors);
    M correctState = UKfilterMath.correctState(Z);
    M correctCov = UKfilterMath.correctCov();
    return correctState;
}
