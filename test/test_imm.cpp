#include <iostream>
#include <catch2/catch.hpp>
#include "ukf.h"
#include "converter.h"
#include "ifilter.h"

using namespace Catch::Benchmark;

template<class M>
struct ConteinerCVCACTxy
{
    Converter<M> converter;
    std::shared_ptr <IFilter<M>> ukfCvSph;
    std::shared_ptr <IFilter<M>> ukfCtSph;
    std::shared_ptr <IFilter<M>> ukfCaSph;

    std::vector<std::shared_ptr<IFilter<M>>> filters;

    void initConteiner(const M& detectionPoint)
    {
        typedef Eigen::SparseMatrix<double> SpMat;
        typedef Eigen::Triplet<double> T;

        M X0;
        M procNoise;
        M measNoise;
        ParamSigmaPoints paramSigmaPoints;

        SpMat Hp(3,6);
        std::vector<T> tripletList;
        tripletList.reserve(3);

        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 2, 1.0));
        tripletList.push_back(T(2, 4, 1.0));
        Hp.setFromTriplets(tripletList.begin(), tripletList.end());
        
        X0 = Hp.transpose() * Utils<M>::sph2CartMeas(detectionPoint);

//------------------------------------------------------
        double process_var;
        double sko_range  = 100.0;
        double sko_Az = 0.1/3.0;
        double sko_Um = 0.1/3.0;
//------------------------------------------------------
        

        process_var = 0.00001;
        procNoise.resize(3, 3);
        procNoise <<    process_var, 0.0, 0.0,
                        0.0, process_var, 0.0,
                        0.0, 0.0, process_var;

        measNoise.resize(3,3);
        measNoise <<   pow(sko_range,2),        0.0,                    0.0,
                            0.0,            pow(sko_Az,2),              0.0,
                            0.0,                0.0,            pow(sko_Um,2);

        paramSigmaPoints.alpha = 1e-3;
        paramSigmaPoints.beta = 2.0;

        paramSigmaPoints.kappa = 3.0 - X0.rows();
        ukfCvSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>>(X0, procNoise, measNoise, paramSigmaPoints);


        procNoise.resize (4, 4);
        procNoise <<    10.0, 0.0, 0.0, 0.0,
                         0.0, 10.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 0.0, 0.0, 1e-7;

        paramSigmaPoints.kappa = -4.0;
        ukfCtSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstTurnXY, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCtXy)}](X0), procNoise, measNoise, paramSigmaPoints);

        process_var = 10.0;
        procNoise.resize(3, 3);
        procNoise <<    process_var, 0.0, 0.0,
                        0.0, process_var, 0.0,
                        0.0, 0.0, process_var;

        paramSigmaPoints.kappa = -6.0;
        ukfCaSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](X0), procNoise, measNoise, paramSigmaPoints);

        filters.push_back(ukfCvSph);
        filters.push_back(ukfCtSph);
        filters.push_back(ukfCaSph);
    }
};


template <class M>
struct StateCov
{
    M state;
    M covariance;
};

template <class M, template <typename> class ConteinerType>
struct Imm
{
    ConteinerType<M> conteiner;
    Converter<M> converter;
    
    M mu_i; // Вероятности режима i
    M p_ij;   // переходная вероятность режима из i в j
    M mu_ij; //смешенная вероятность
    M cj;
    
    Imm(const M& modeProbability, const M& transmitProbability, const ConteinerType<M>& conteiner):
    conteiner(conteiner), mu_i(modeProbability), p_ij(transmitProbability) {}

    M step(const M& Z ,double dt)
    {
        PRINTM (mu_i);
        mu_ij = computeMixingProbability(p_ij, mu_i);
        PRINTM(mu_ij);
        std::pair<std::vector<M>,std::vector<M>> stateCovInit = InitMixingStateAndCovariance (mu_ij);

        filterStep(Z, stateCovInit,dt);
        updateModeProbability(Z, cj);
        return combinationModelCondition();

    }
    M step(double dt)
    {
        mu_ij = computeMixingProbability(p_ij, mu_i);
        std::pair<std::vector<M>, std::vector<M>> stateCovInit = InitMixingStateAndCovariance(mu_ij);
        for (size_t i = 0; i < conteiner.filters.size(); ++i)
        {   
            conteiner.filters[i]->correctStruct.X = stateCovInit.first[i]; // допустим тут будет своп
            conteiner.filters[i]->correctStruct.P = stateCovInit.second[i];
    
            conteiner.filters[i]->correctStruct.X = conteiner.filters[i]->predict(dt);
            conteiner.filters[i]->correctStruct.P = conteiner.filters[i]->predictStruct.Pe;

            mu_i(0, i) = cj(0, i);
        }

        return combinationModelCondition();
    }

    M computeMixingProbability(const M &p_ij, const M &mu_i)
    {
        cj.resize(1, mu_i.cols());

        for (long int j = 0; j < mu_i.cols(); j++)
        {
            double c = 0.;
            for (long int i = 0; i < p_ij.cols(); i++)
            {
                c +=(p_ij(i, j) * mu_i(0, i));
            }
            cj(0, j) = c;
        }

        M mix(p_ij.rows(), p_ij.cols());

        for (long int j = 0; j < p_ij.cols(); j++)
        {
            for (long int i = 0; i < p_ij.rows(); i++)
            {
                mix(i, j) = p_ij(i, j) * mu_i(0, i) / cj(0, j);
            }
        }
        return mix;
    }


       std::pair<std::vector<M>, std::vector<M>> InitMixingStateAndCovariance(const M &mu_ij)
    {
        std::vector<M> statesOfFilters;
        std::vector<M> covarianceOfFilters;

        for (size_t j = 0; j < conteiner.filters.size(); ++j)
        {
            M mixedState = M::Zero(conteiner.filters[j]->correctStruct.X.rows(), conteiner.filters[j]->correctStruct.X.cols());
            M mixedCovariance = M::Zero(conteiner.filters[j]->correctStruct.P.rows(), conteiner.filters[j]->correctStruct.P.cols());

            for (size_t i = 0; i < conteiner.filters.size(); ++i)
            {
                M convertedState = converter.m[{conteiner.filters[i]->getModelType(),conteiner.filters[j]->getModelType()}](conteiner.filters[i]->correctStruct.X);
                mixedState += mu_ij(i, j) * convertedState;
            }
                statesOfFilters.push_back(mixedState);

            for (size_t i = 0; i < conteiner.filters.size(); ++i)
            {
                M convertedState = converter.m[{conteiner.filters[i]->getModelType(),conteiner.filters[j]->getModelType()}](conteiner.filters[i]->correctStruct.X);
                M dX = convertedState - mixedState;
                M convertedCovariance = converter.m[{conteiner.filters[i]->getModelType(),conteiner.filters[j]->getModelType()}](conteiner.filters[i]->correctStruct.P);
                mixedCovariance += mu_ij(i, j) * (convertedCovariance + dX * dX.transpose());
            }
                        
            covarianceOfFilters.push_back(mixedCovariance);
        }

        return std::make_pair(statesOfFilters, covarianceOfFilters);
    }

    void filterStep(const M &Z, std::pair<std::vector<M>, std::vector<M>> &stateCov, double dt)
    {
        for (size_t i = 0; i < conteiner.filters.size(); ++i)
        {

            conteiner.filters[i]->correctStruct.X = stateCov.first[i];
            conteiner.filters[i]->correctStruct.P = stateCov.second[i];

            conteiner.filters[i]->predict(dt);
            conteiner.filters[i]->correct(Z);
        }
    }

    double likelihoodFunction(const M &Z, const M &Ze, const M &Se)
    {
        double n = Z.rows();
        M v = Z - Ze;
        long double power = -0.5 * (v.transpose() * Se.inverse() * v)(0, 0);
        long double probability = std::pow((1 / (2 * M_PI)), n / 2.0) / std::sqrt(Se.determinant()) * std::exp(power);

        return probability;
    }

    void updateModeProbability(const M &Z, const M &cj)
    {
        long double c = 0.0;
        for (size_t i = 0; i < conteiner.filters.size(); ++i)
        {
            mu_i(0, i) = likelihoodFunction(Z, conteiner.filters[i]->predictStruct.Ze, conteiner.filters[i]->predictStruct.Se) * cj(0, i);
            c += mu_i(0, i);
        }

        for (size_t i = 0; i < conteiner.filters.size(); ++i)
        {
            mu_i(0, i) /= c;
        }

    }

    M combinationModelCondition(/* Флаг означающий модель вывода состояния*/) // сейчас возвращаяется модель CV
    {
        M X = M::Zero(conteiner.filters[0]->correctStruct.X.rows(), conteiner.filters[0]->correctStruct.X.cols());
        M dx = M::Zero(conteiner.filters[0]->correctStruct.X.rows(), conteiner.filters[0]->correctStruct.X.cols());
        M P = M::Zero(conteiner.filters[0]->correctStruct.P.rows(), conteiner.filters[0]->correctStruct.P.cols());

        for (size_t i = 0; i < conteiner.filters.size(); ++i)
        {
            M convertedState = converter.m[{conteiner.filters[i]->getModelType(), typeid(converter.modelCv)}](conteiner.filters[i]->correctStruct.X);
            X += mu_i(0, i) * convertedState;
        }

        for (size_t i = 0; i < conteiner.filters.size(); ++i)
        {
            dx = converter.m[{conteiner.filters[i]->getModelType(), typeid(converter.modelCv)}](conteiner.filters[i]->correctStruct.X) - X;
            M convertedCovariance = converter.m[{conteiner.filters[i]->getModelType(), typeid(converter.modelCv)}](conteiner.filters[i]->correctStruct.P);
            P += mu_i(0, i) * (convertedCovariance + dx * dx.transpose());
        }
        return X;
    }
};

TEST_CASE("test_imm")
{
    Eigen::MatrixXd X(6, 1);
    Eigen::MatrixXd Z0(3, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Z1(3, 1);
    Eigen::MatrixXd Q(3, 3);
    Eigen::MatrixXd R(3, 3);

    
    Z0 << 2.43146250e+04, 6.30715510e+01, 2.40279867e+01; 
    Z << 2.44146250e+04, 6.36715510e+01, 2.41279867e+01;      
    Z1 << 2.45146250e+04, 6.37715510e+01, 2.42279867e+01;                
 
    Q << 10.0,0.0,0.0,
        0.0,10.0,0.0,
        0.0,0.0,10.0;
        
    R << 10000.0, 0.0, 0.0,
        0.0, pow((0.1/3),2), 0.0,
        0.0, 0.0, pow((0.1/3),2);
    
    ParamSigmaPoints p;
    p.alpha = 1e-3;
    p.beta = 2;
    p.kappa = -3;


    Eigen::MatrixXd mui(1,3);
    Eigen::MatrixXd Pij(3,3);
    mui <<1.0/3.0,1.0/3.0,1.0/3.0;

    Pij <<   0.97, 0.015, 0.015,
            0.015,  0.97, 0.015,
            0.015, 0.015,  0.97;


    ConteinerCVCACTxy<Eigen::MatrixXd> conteiner;
    conteiner.initConteiner(Z0);
    Imm<Eigen::MatrixXd, ConteinerCVCACTxy> imm(mui, Pij, conteiner);

    double dt = 1.0;
    PRINTM(imm.step(Z,dt));
    // PRINTM(imm.step(dt));
    PRINTM(imm.step(Z1,dt));
    BENCHMARK("imm"){

    };
}
