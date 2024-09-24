#pragma once

#include <vector>
#include <iostream>
#include "utils.h"
#include "models.h"
#include "ukfMath.h"
#include "ifilter.h"

template <class M,
          class TypeInitialization,
          class TypeStateModelFunc,
          class TypeMeasFunc,
          class TypeControlFunc>
struct UnscentedKalmanFilter
    : public IFilter<M>
{
    std::pair<M, M> predict(double dt) override final
    {
        UKfilterMath.compute_weights(paramsSigmaPoints, lamda, c, correctInfo.X.rows(), Wc, Wm);
        UKfilterMath.compute_sigma_points(correctInfo.X, correctInfo.P, lamda, sigmaVectors, U);
        extrapolatedStateSigmaVectors = stateFunc(sigmaVectors, dt);
        UKfilterMath.doExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictInfo.Xe, Wm);
        G = controlFunc(dt);
        UKfilterMath.doCovMatExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictInfo.Xe, G, Q, predictInfo.Pe, Wc);
        correctInfo.X = predictInfo.Xe;  // Записываю в качестве скореектированных предсказанные значения. под вопросом! Нужно в track.step(dt)
        correctInfo.P = predictInfo.Pe;

        return std::make_pair(predictInfo.Xe, predictInfo.Pe);
    }

    std::pair<M, M> correct(const M &Z) override final
    {
        extrapolatedMeasSigmaVectors = measFunc(extrapolatedStateSigmaVectors, Z);
        UKfilterMath.doExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictInfo.Ze, Wm);
        UKfilterMath.doCovMatExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictInfo.Ze, predictInfo.Pzz, Wc);
        UKfilterMath.doCovMatInnovation(predictInfo.Pzz, R, predictInfo.Se);
        UKfilterMath.calcGainFilter(extrapolatedStateSigmaVectors, predictInfo.Xe, extrapolatedMeasSigmaVectors, predictInfo.Ze, predictInfo.Se, predictInfo.Pxz, predictInfo.K, Wc);
        UKfilterMath.correctState(predictInfo.Xe, Z, predictInfo.Ze, predictInfo.K, correctInfo.X);
        UKfilterMath.correctCov(predictInfo.Pe, predictInfo.K, predictInfo.Se, correctInfo.P);
        return std::make_pair(correctInfo.X, correctInfo.P);
    }

    double likelihood(const M &Z) override final
    {
        v = Z - predictInfo.Ze;
        long double power = -0.5 * (v.transpose() * predictInfo.Se.inverse() * v)(0, 0);
        long double probability = std::pow((1 / (2 * M_PI)), Z.rows() / 2.0) / std::sqrt(predictInfo.Se.determinant()) * std::exp(power);
        return probability;
    }

    double distance(const M &Z) override final
    {

        v = Z - predictInfo.Ze;
        double distance = (v.transpose() * predictInfo.Se.inverse() * v)(0,0);
        return distance;
    }

    std::type_index getModelType() const override
    {
        return std::type_index(typeid(stateFunc));
    }

    Correct<M> getCorrectInfo() override final {

        return correctInfo;
    }

    Predict<M> getPredictInfo() override final {

        return predictInfo;
    }
    void setCorrectInfo(const M& X, const M& P) override final
    {
        correctInfo.X = X;
        correctInfo.P = P;
    }

    UnscentedKalmanFilter()
    {
        predictInfo.Xe.resize(stateFunc.getSize(), 1);
        predictInfo.Pe.resize(stateFunc.getSize(), stateFunc.getSize());
        predictInfo.Ze.resize(measFunc.getSize(), 1);
        predictInfo.Pzz.resize(measFunc.getSize(), measFunc.getSize());
        predictInfo.Pxz.resize(stateFunc.getSize(), measFunc.getSize());
        predictInfo.Se.resize(measFunc.getSize(), measFunc.getSize());
        predictInfo.K.resize(stateFunc.getSize(), measFunc.getSize());
        correctInfo.X.resize(stateFunc.getSize(), 1);
        correctInfo.P.resize(stateFunc.getSize(), stateFunc.getSize());
        Q.resize(measFunc.getSize(), measFunc.getSize());
        R.resize(measFunc.getSize(), measFunc.getSize());
        G.resize(controlFunc.getSize().first, controlFunc.getSize().second);
        v.resize(measFunc.getSize(), 1);
        sigmaVectors.resize(stateFunc.getSize(), 2 * stateFunc.getSize() + 1);
        Wc.reserve(sigmaVectors.cols());
        Wm.reserve(sigmaVectors.cols());
        U.resize(correctInfo.P.rows(), correctInfo.P.cols());
        extrapolatedStateSigmaVectors.resize(sigmaVectors.rows(), sigmaVectors.cols());
        extrapolatedMeasSigmaVectors.resize(measFunc.getSize(), extrapolatedStateSigmaVectors.cols());
    }

    template <class... TypeArgs>
    void Initialization(TypeArgs... args)
    {
        initializator(*this, args...);
    }

    Predict<M> predictInfo;
    Correct<M> correctInfo;
    ParamSigmaPoints paramsSigmaPoints;
    TypeStateModelFunc stateFunc;
    TypeMeasFunc measFunc;
    M Q;
    M R;

private:
    TypeInitialization initializator;

    TypeControlFunc controlFunc;
    UnscentedKalmanFilterMath<M> UKfilterMath;
    M extrapolatedStateSigmaVectors;
    M extrapolatedMeasSigmaVectors;
    M G;
    M v;
    M sigmaVectors;
    M U;
    double lamda;
    double c;
    std::vector<double> Wc, Wm;
};