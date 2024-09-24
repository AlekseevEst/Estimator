#pragma once
#include <iostream>
#include "converter.h"
#include "ifilter.h"
#include "immMath.h"

template <class M,
          class TypeInitialization> 
struct IMM
    : public IFilter<M>
{
    std::pair<M, M> predict(double dt) override final
    {
        math.computeMixingProbability(p_ij, mu_i, mu_ij, cj);
        math.MixingStateAndCovariance(mu_ij, stateMixed, covarianceMixed, initializator.filters, converter);
        for (size_t i = 0; i < initializator.filters.size(); ++i)
            {
                initializator.filters[i]->setCorrectInfo(stateMixed[i],covarianceMixed[i]);
                initializator.filters[i]->predict(dt);
                
            }
        return math.combinationModelsCondition(mu_i, predictInfo.Xe, predictInfo.Pe, initializator.filters, converter);
    }
      
    std::pair<M, M> correct(const M &Z) override final
    {
        for (size_t i = 0; i < initializator.filters.size(); ++i)
        {
            auto cor = initializator.filters[i]->correct(Z);
        }

        math.updateModeProbability(Z, cj, mu_i, initializator.filters);
        return math.combinationModelsCondition(mu_i, correctInfo.X, correctInfo.P, initializator.filters, converter);
    }

    double likelihood(const M &Z) override final
    { 
        double probability = 0.;
        for (size_t i = 0; i < initializator.filters.size(); ++i)
        {
            probability += mu_i(0, i) * initializator.filters[i]->likelihood(Z);
        }
        return probability;
    }

    double distance(const M &Z) override final
    {
        double totalDistance = 0.0;
        for (size_t i = 0; i < initializator.filters.size(); ++i)
        {
        M v = Z - initializator.filters[i]->getPredictInfo().Ze;
        double mahalonobisDistance = (v.transpose() * initializator.filters[i]->getPredictInfo().Se.inverse() * v)(0,0);
        totalDistance += mahalonobisDistance * mu_i(0,i);
        }
        return totalDistance;
    }

    std::type_index getModelType() const override
    {
        return std::type_index(typeid(initializator));
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

    IMM()
    {
        mu_i.resize(1, initializator.filters.size());
        p_ij.resize(initializator.filters.size(),initializator.filters.size());
        mu_ij.resize(p_ij.rows(),p_ij.cols());
        cj.resize(1, mu_i.cols());
        stateMixed.resize(initializator.filters.size());
        covarianceMixed.resize(initializator.filters.size());

        correctInfo.X.resize(initializator.filterCV->correctInfo.X.rows(), initializator.filterCV->correctInfo.X.cols());
        correctInfo.P.resize(correctInfo.X.rows(), correctInfo.X.rows());
        predictInfo.Xe.resize(initializator.filterCV->correctInfo.X.rows(), initializator.filterCV->correctInfo.X.cols());
        predictInfo.Pe.resize(correctInfo.X.rows(), correctInfo.X.rows());


        for (size_t i = 0; i < initializator.filters.size(); ++i)
        {
            stateMixed[i].resize(initializator.filters[i]->getCorrectInfo().X.rows(),initializator.filters[i]->getCorrectInfo().X.cols());
            covarianceMixed[i].resize(stateMixed[i].rows(),stateMixed[i].rows());
        }    

    }

    template <class... TypeArgs>
    void Initialization(TypeArgs... args)
    {
        initializator(*this, args...);
    }
    M mu_ij; // смешенная вероятность
    M mu_i; // Вероятности режима i
    M p_ij; // переходная вероятность режима из i в j
    M cj;
    Predict<M> predictInfo;
    Correct<M> correctInfo;

    std::vector<M> stateMixed;
    std::vector<M> covarianceMixed;
    TypeInitialization initializator;
private:

    ImmMath<M> math;
    Converter<M> converter;

   
   


};