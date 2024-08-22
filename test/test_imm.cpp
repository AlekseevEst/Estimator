#include <iostream>
#include <catch2/catch.hpp>
#include "ukf.h"
#include "converter.h"
#include "ifilter.h"

using namespace Catch::Benchmark;

template<class M>
struct ConteinerCVCACT
{
    Converter<M> converter;
    std::shared_ptr <IFilter<M>> ukfCvSph;
    std::shared_ptr <IFilter<M>> ukfCtSph;
    std::shared_ptr <IFilter<M>> ukfCaSph;

    ConteinerCVCACT(const M& X, const M& procNoise, const M& measNoise, ParamSigmaPoints paramSigmaPoints)
    {
               M ProcNoiseCT(4, 4);
        ProcNoiseCT <<   10.0, 0.0, 0.0, 0.0,
                         0.0, 10.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 0.0, 0.0, 1e-7;

        M ProcNoiseCV(3, 3);
        ProcNoiseCV <<   0.00001, 0.0, 0.0,
                         0.0, 0.00001, 0.0,
                         0.0, 0.0, 0.00001;


        // здесь выделяем память. Пока что все вместе: и выделение и инициализация фильтров
        ukfCvSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>>(X, ProcNoiseCV, measNoise, paramSigmaPoints);
        paramSigmaPoints.kappa = 3 - 7; //параметр ансцентного преобразования задал вручную для CT. 
        ukfCtSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstTurnXZ, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCt)}](X), ProcNoiseCT, measNoise, paramSigmaPoints);
        paramSigmaPoints.kappa = 3 - 9; //параметр ансцентного преобразования задал вручную для CA. 
        ukfCaSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](X), procNoise, measNoise, paramSigmaPoints);

    }
    // initConteinerCVCACT()
    // {
        
    //     // инициализиция фильтров своим состоянием. ????????
    // }
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

    // StateCov<M> stateCovCorrect;
    // StateCov<M> statCovBuff;
    // StateCov<M> stateCovPredict;
    
    M mu_i; // Вероятности режима i
    M p_ij;   // переходная вероятность режима из i в j
    M mu_ij; //смешенная вероятность
    M cj;
    
    Imm(const M& modeProbability, const M& transmitProbability, const M& X, const M& procNoise, const M& measNoise, ParamSigmaPoints paramSigmaPoints ):
    conteiner(X,procNoise,measNoise,paramSigmaPoints)
    {
        mu_i = modeProbability;
        p_ij = transmitProbability;
    }
    // void immInit()

    M step(const M& Z ,double dt)
    {
        mu_ij = computeMixingProbability(p_ij, mu_i);
        PRINTM(mu_ij);
        std::pair<std::vector<M>,std::vector<M>> stateCovInit = InitMixingStateAndCovariance (mu_ij);

        filterStep(Z, stateCovInit,dt);
        updateModeProbability(Z, cj);
        M X = combinationModelCondition();

        return X;
    }

     M step (double dt)
    {
        mu_ij = computeMixingProbability(p_ij, mu_i);
        std::pair<std::vector<M>, std::vector<M>> stateCovInit = InitMixingStateAndCovariance(mu_ij);

        conteiner.ukfCvSph->correctStruct.X = stateCovInit.first[0];
        conteiner.ukfCvSph->correctStruct.P = stateCovInit.second[0];
        conteiner.ukfCvSph->correctStruct.X = conteiner.ukfCvSph->predict(dt);
        conteiner.ukfCvSph->correctStruct.P = conteiner.ukfCvSph->predictStruct.Pe;

        conteiner.ukfCtSph->correctStruct.X = stateCovInit.first[1];
        conteiner.ukfCtSph->correctStruct.P = stateCovInit.second[1];
        conteiner.ukfCtSph->correctStruct.X = conteiner.ukfCtSph->predict(dt);
        conteiner.ukfCtSph->correctStruct.P = conteiner.ukfCtSph->predictStruct.Pe;

        conteiner.ukfCaSph->correctStruct.X = stateCovInit.first[2];
        conteiner.ukfCaSph->correctStruct.P = stateCovInit.second[2];
        conteiner.ukfCaSph->correctStruct.X = conteiner.ukfCaSph->predict(dt);
        conteiner.ukfCaSph->correctStruct.P = conteiner.ukfCaSph->predictStruct.Pe;

        mu_i(0,0) = cj(0,0);
        mu_i(0,1) = cj(0,1);
        mu_i(0,2) = cj(0,2);

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
                c = c + (p_ij(i, j) * mu_i(0, i));
            }
            cj(0, j) = c;
        }

        M mix;
        mix.resize(p_ij.rows(), p_ij.cols());

        for (long int j = 0; j < p_ij.cols(); j++)
        {
            for (long int i = 0; i < p_ij.rows(); i++)
            {
                mix(i, j) = p_ij(i, j) * mu_i(0, i) / cj(0, j);
            }
        }
        return mix;
    }

    std::pair<std::vector<M>,std::vector<M>> InitMixingStateAndCovariance(const M& mu_ij) // как сделать универсально?
    {
        std::pair<std::vector<M>,std::vector<M>> initMixingStateAndCovariance;
        std::vector<M> statesOfFilters;
        std::vector<M> covarianceOfFilters;

        // PRINTM(mu_ij);
        // PRINTM(conteiner.ukfCvSph->correctStruct.P);
        // PRINTM(conteiner.ukfCtSph->correctStruct.P);
        // PRINTM(conteiner.ukfCaSph->correctStruct.P);


        M X0Cv = mu_ij(0,0) * conteiner.ukfCvSph->correctStruct.X + 
                 mu_ij(1,0) * converter.m[{typeid(converter.modelCt), typeid(converter.modelCv)}](conteiner.ukfCtSph->correctStruct.X) + 
                 mu_ij(2,0) * converter.m[{typeid(converter.modelCa), typeid(converter.modelCv)}](conteiner.ukfCaSph->correctStruct.X);

        statesOfFilters.push_back(X0Cv);

        M X0Ct = mu_ij(0,1) * converter.m[{typeid(converter.modelCv), typeid(converter.modelCt)}](conteiner.ukfCvSph->correctStruct.X) + 
                 mu_ij(1,1) * conteiner.ukfCtSph->correctStruct.X + 
                 mu_ij(2,1) * converter.m[{typeid(converter.modelCa), typeid(converter.modelCt)}](conteiner.ukfCaSph->correctStruct.X);  
  
        statesOfFilters.push_back(X0Ct);

        M X0Ca = mu_ij(0,2) * converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](conteiner.ukfCvSph->correctStruct.X) + 
                 mu_ij(1,2) * converter.m[{typeid(converter.modelCt), typeid(converter.modelCa)}](conteiner.ukfCtSph->correctStruct.X) + 
                 mu_ij(2,2) * conteiner.ukfCaSph->correctStruct.X;

        statesOfFilters.push_back(X0Ca);

        M dXcvX0cv = conteiner.ukfCvSph->correctStruct.X - X0Cv;
        M dXctX0cv = converter.m[{typeid(converter.modelCt), typeid(converter.modelCv)}](conteiner.ukfCtSph->correctStruct.X) - X0Cv;
        M dXcaX0cv = converter.m[{typeid(converter.modelCa), typeid(converter.modelCv)}](conteiner.ukfCaSph->correctStruct.X) - X0Cv;


        M dXcvX0ct = converter.m[{typeid(converter.modelCv), typeid(converter.modelCt)}](conteiner.ukfCvSph->correctStruct.X) - X0Ct;
        M dXctX0ct = conteiner.ukfCtSph->correctStruct.X - X0Ct;
        M dXcaX0ct = converter.m[{typeid(converter.modelCa), typeid(converter.modelCt)}](conteiner.ukfCaSph->correctStruct.X) - X0Ct;


        M dXcvX0ca = converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](conteiner.ukfCvSph->correctStruct.X) - X0Ca;                          
        M dXctX0ca = converter.m[{typeid(converter.modelCt), typeid(converter.modelCa)}](conteiner.ukfCtSph->correctStruct.X) - X0Ca;        
        M dXcaX0ca = conteiner.ukfCaSph->correctStruct.X - X0Ca;  


        M P0Cv = mu_ij(0,0) * (conteiner.ukfCvSph->correctStruct.P + dXcvX0cv * dXcvX0cv.transpose()) +
                 mu_ij(1,0) * (converter.m[{typeid(converter.modelCt), typeid(converter.modelCv)}](conteiner.ukfCtSph->correctStruct.P) + dXctX0cv * dXctX0cv.transpose()) +
                 mu_ij(2,0) * (converter.m[{typeid(converter.modelCa), typeid(converter.modelCv)}](conteiner.ukfCaSph->correctStruct.P) + dXcaX0cv * dXcaX0cv.transpose());
        covarianceOfFilters.push_back(P0Cv);

        M P0Ct = mu_ij(0,1) * (converter.m[{typeid(converter.modelCv), typeid(converter.modelCt)}](conteiner.ukfCvSph->correctStruct.P) + dXcvX0ct * dXcvX0ct.transpose()) +
                 mu_ij(1,1) * (conteiner.ukfCtSph->correctStruct.P + dXctX0ct * dXctX0ct.transpose()) +
                 mu_ij(2,1) * (converter.m[{typeid(converter.modelCa), typeid(converter.modelCt)}](conteiner.ukfCaSph->correctStruct.P) + dXcaX0ct * dXcaX0ct.transpose());
        covarianceOfFilters.push_back(P0Ct);
                 
        M P0Ca = mu_ij(0,2) * (converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](conteiner.ukfCvSph->correctStruct.P) + dXcvX0ca * dXcvX0ca.transpose()) +
                 mu_ij(1,2) * (converter.m[{typeid(converter.modelCt), typeid(converter.modelCa)}](conteiner.ukfCtSph->correctStruct.P) + dXctX0ca * dXctX0ca.transpose()) +
                 mu_ij(2,2) * (conteiner.ukfCaSph->correctStruct.P + dXcaX0ca * dXcaX0ca.transpose());    
        covarianceOfFilters.push_back(P0Ca);
        
        initMixingStateAndCovariance = std::make_pair(statesOfFilters,covarianceOfFilters);

        return initMixingStateAndCovariance;
    }

    void filterStep(const M& Z, std::pair<std::vector<M>,std::vector<M>>& stateCov, double dt)
    {   
        conteiner.ukfCvSph->correctStruct.X = stateCov.first[0];
        conteiner.ukfCvSph->correctStruct.P = stateCov.second[0];

        conteiner.ukfCvSph->predict(dt);
        conteiner.ukfCvSph->correct(Z);

        conteiner.ukfCtSph->correctStruct.X = stateCov.first[1];
        conteiner.ukfCtSph->correctStruct.P = stateCov.second[1];

        conteiner.ukfCtSph->predict(dt);
        conteiner.ukfCtSph->correct(Z);

        conteiner.ukfCaSph->correctStruct.X = stateCov.first[2];
        conteiner.ukfCaSph->correctStruct.P = stateCov.second[2];

        conteiner.ukfCaSph->predict(dt);
        conteiner.ukfCaSph->correct(Z);   
    }

    double likelihoodFunction(const M &Z, const M &Ze, const M &Se)
    {
        double n = Z.rows();
        M v = Z - Ze;
        long double power = -0.5 * (v.transpose() * Se.inverse() * v)(0, 0);
        long double probability = std::pow((1 / (2 * M_PI)), n / 2.0) / std::sqrt(Se.determinant()) * std::exp(power);

        return probability;
    }

    void updateModeProbability(const M& Z, const M& cj)
    {
        mu_i(0,0) = likelihoodFunction(Z, conteiner.ukfCvSph->predictStruct.Ze, conteiner.ukfCvSph->predictStruct.Se) * cj(0,0);
        mu_i(0,1) = likelihoodFunction(Z, conteiner.ukfCtSph->predictStruct.Ze, conteiner.ukfCtSph->predictStruct.Se) * cj(0,1);
        mu_i(0,2) = likelihoodFunction(Z, conteiner.ukfCaSph->predictStruct.Ze, conteiner.ukfCaSph->predictStruct.Se) * cj(0,2);
        long double c = mu_i(0,0) + mu_i(0,1) + mu_i(0,2);

        mu_i(0,0) = mu_i(0,0)/c;
        mu_i(0,1) = mu_i(0,1)/c;
        mu_i(0,2) = mu_i(0,2)/c;
    }

    M combinationModelCondition()
    {
       M X = mu_i(0,0)* conteiner.ukfCvSph->correctStruct.X +
             mu_i(0,1)* converter.m[{typeid(converter.modelCt), typeid(converter.modelCv)}] (conteiner.ukfCtSph->correctStruct.X) +
             mu_i(0,2)* converter.m[{typeid(converter.modelCa), typeid (converter.modelCv)}] (conteiner.ukfCaSph->correctStruct.X);

        M dxCv = conteiner.ukfCvSph->correctStruct.X - X;
        M dxCt = converter.m[{typeid(converter.modelCt), typeid(converter.modelCv)}] (conteiner.ukfCtSph->correctStruct.X) - X;
        M dxCa = converter.m[{typeid(converter.modelCa), typeid(converter.modelCv)}] (conteiner.ukfCaSph->correctStruct.X) - X;

       M P =  mu_i(0,0) * (conteiner.ukfCvSph->correctStruct.P + dxCv * dxCv.transpose()) +
              mu_i(0,1) * (converter.m[{typeid(converter.modelCt), typeid(converter.modelCv)}](conteiner.ukfCtSph->correctStruct.P) + dxCt * dxCt.transpose())+
              mu_i(0,2) * (converter.m[{typeid(converter.modelCa), typeid(converter.modelCv)}](conteiner.ukfCaSph->correctStruct.P) + dxCa * dxCa.transpose());

        return X;
    }
};

TEST_CASE("test_imm")
{
    Eigen::MatrixXd X(6, 1);
    Eigen::MatrixXd Z(3, 1);
    Eigen::MatrixXd Q(3, 3);
    Eigen::MatrixXd R(3, 3);

    X << 9978.99978886, 0.0, 20033.44840212, 0.0, 10033.56458559, 0.0;
    Z << 2.44146250e+04, 6.36715510e+01, 2.41279867e+01;                 
 
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



    Imm<Eigen::MatrixXd, ConteinerCVCACT> imm(mui, Pij, X, Q, R, p);

    double dt = 1.0;
    PRINTM(imm.step(Z,dt));
    PRINTM(imm.step(dt));
    BENCHMARK("imm"){

    };
}
