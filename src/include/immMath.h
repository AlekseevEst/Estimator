#pragma once
#include "vector"
#include "ifilter.h"
template <class M> 
struct ImmMath{

    void computeMixingProbability(const M &p_ij, const M &mu_i, M &mu_ij, M &cj)
    {
        cj = mu_i * p_ij;
        for (long int j = 0; j < p_ij.cols(); j++)
        {
            for (long int i = 0; i < p_ij.rows(); i++)
            {
                mu_ij(i, j) = p_ij(i, j) * mu_i(0, i) / cj(0, j);
            }
        }
    }

    void MixingStateAndCovariance(const M &mu_ij, std::vector<M> &stateMixed, std::vector<M> &covarianceMixed, const std::vector<std::shared_ptr<IFilter<M>>> &filters, Converter<M> &converter)
    {
        for (size_t j = 0; j < filters.size(); ++j)
        {
            stateMixed[j].setZero();
            covarianceMixed[j].setZero();
        }

        for (size_t j = 0; j < filters.size(); ++j)
        {
            for (size_t i = 0; i < filters.size(); ++i)

            {
                M convertedState = converter.m[{filters[i]->getModelType(), filters[j]->getModelType()}](filters[i]->getCorrectInfo().X);
                stateMixed[j] += mu_ij(i, j) * convertedState;
            }
            for (size_t i = 0; i < filters.size(); ++i)
            {                
                M convertedState = converter.m[{filters[i]->getModelType(),filters[j]->getModelType()}](filters[i]->getCorrectInfo().X);
                M dX = convertedState - stateMixed[j];
                M convertedCovariance = converter.m[{filters[i]->getModelType(), filters[j]->getModelType()}](filters[i]->getCorrectInfo().P);
                covarianceMixed[j] += mu_ij(i, j) * (convertedCovariance + dX * dX.transpose());
            }
        }
    }

    void updateModeProbability(const M &Z, const M &cj, M& mu_i, std::vector<std::shared_ptr<IFilter<M>>>& filters)
    {
        long double c = 0.0;

        for (size_t i = 0; i < filters.size(); ++i)
        {
            mu_i(0, i) = filters[i]->likelihood(Z) * cj(0, i);
            c += mu_i(0, i);
        }

        for (size_t i = 0; i < filters.size(); ++i)
        {
            mu_i(0, i) /= c;
        }

    }

        std::pair<M,M> combinationModelsCondition(M& mu_i, M& X, M& P, std::vector<std::shared_ptr<IFilter<M>>>& filters, Converter<M>& converter) 
    {
        X.setZero();
        P.setZero();
        for (size_t i = 0; i < filters.size(); ++i)
        {
            M convertedState = converter.m[{filters[i]->getModelType(), filters[0]->getModelType()}](filters[i]->getCorrectInfo().X);
            X += mu_i(0, i) * convertedState;
        }

        for (size_t i = 0; i < filters.size(); ++i)
        {
            M dx = converter.m[{filters[i]->getModelType(), filters[0]->getModelType()}](filters[i]->getCorrectInfo().X) - X;
            M convertedCovariance = converter.m[{filters[i]->getModelType(), filters[0]->getModelType()}](filters[i]->getCorrectInfo().P);
            P += mu_i(0, i) * (convertedCovariance + dx * dx.transpose());
        }

        return std::make_pair(X,P);
    }


private:

};