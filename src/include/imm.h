#pragma once
#include <iostream>
#include "ukf.h"
#include "converter.h"
#include "conteiner.h"
#include "ifilter.h"

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
    {   std::cout<< "CV" << std::endl<< std::endl;
        // std::cout<< "На это будем метять coorectStruct.X" << std::endl<< std::endl;
        // PRINTM(stateCov.first[0]);
        // std::cout<< "Было coorectStruct.X" << std::endl<< std::endl;
        // PRINTM(conteiner.ukfCvSph->correctStruct.X);
        conteiner.ukfCvSph->correctStruct.X = stateCov.first[0];
        // std::cout<< "стало coorectStruct.X" << std::endl<< std::endl;
        // PRINTM(conteiner.ukfCvSph->correctStruct.X);
        
        // std::cout<< "На это будем метять coorectStruct.Р" << std::endl<< std::endl;
        // PRINTM(stateCov.second[0]);
        // std::cout<< "Было coorectStruct.P" << std::endl<< std::endl;
        conteiner.ukfCvSph->correctStruct.P = stateCov.second[0];
        // std::cout<< "стало coorectStruct.Р" << std::endl<< std::endl;

        // std::cout<< "После Predict_CV" << std::endl<< std::endl;

        conteiner.ukfCvSph->predict(dt);

        // PRINTM(conteiner.ukfCvSph->correctStruct.X);
        // PRINTM(conteiner.ukfCvSph->correctStruct.P);

        std::cout<< "После Correct_CV" << std::endl<< std::endl;
        conteiner.ukfCvSph->correct(Z);
        // PRINTM(conteiner.ukfCvSph->correctStruct.X);
        // PRINTM(conteiner.ukfCvSph->correctStruct.P);


        std::cout<< "CT" << std::endl<< std::endl;
        // std::cout<< "Было coorectStruct.X и Р" << std::endl<< std::endl;
        // PRINTM(conteiner.ukfCtSph->correctStruct.X);
        // PRINTM(conteiner.ukfCtSph->correctStruct.P);
        // std::cout<< "На это будем метять coorectStruct.X и Р" << std::endl<< std::endl;
        // PRINTM(stateCov.first[1]);
        // PRINTM(stateCov.second[1]);
        conteiner.ukfCtSph->correctStruct.X = stateCov.first[1];
        conteiner.ukfCtSph->correctStruct.P = stateCov.second[1];

        // std::cout<< "Стало coorectStruct.X и Р" << std::endl<< std::endl;
        // PRINTM(conteiner.ukfCtSph->correctStruct.X);
        // PRINTM(conteiner.ukfCtSph->correctStruct.P);

        // std::cout<< "Полсе Predict_CT" << std::endl<< std::endl;
        conteiner.ukfCtSph->predict(dt);
        // PRINTM(conteiner.ukfCtSph->correctStruct.X);
        // PRINTM(conteiner.ukfCtSph->correctStruct.P);

        std::cout<< "После Correct_CT" << std::endl<< std::endl;
        conteiner.ukfCtSph->correct(Z);
        // PRINTM(conteiner.ukfCtSph->correctStruct.X);
        // PRINTM(conteiner.ukfCtSph->correctStruct.P);

        std::cout<< "CA" << std::endl<< std::endl;
        // std::cout<< "Было coorectStruct.X и Р" << std::endl<< std::endl;
        // PRINTM(conteiner.ukfCaSph->correctStruct.X);
        // PRINTM(conteiner.ukfCaSph->correctStruct.P);
        // std::cout<< "На это будем метять coorectStruct.X и Р" << std::endl<< std::endl;
        // PRINTM(stateCov.first[2]);
        // PRINTM(stateCov.second[2]);
        conteiner.ukfCaSph->correctStruct.X = stateCov.first[2];
        conteiner.ukfCaSph->correctStruct.P = stateCov.second[2];
        // std::cout<< "Стало coorectStruct.X и Р" << std::endl<< std::endl;

        // PRINTM(conteiner.ukfCaSph->correctStruct.X);
        // PRINTM(conteiner.ukfCaSph->correctStruct.P);

        // std::cout<< "После Predict_CA" << std::endl<< std::endl;
        conteiner.ukfCaSph->predict(dt);

        // PRINTM(conteiner.ukfCaSph->correctStruct.X);
        // PRINTM(conteiner.ukfCaSph->correctStruct.P);

        std::cout<< "После Correct_CA" << std::endl<< std::endl;
        conteiner.ukfCaSph->correct(Z);
        // PRINTM(conteiner.ukfCaSph->correctStruct.X);
        // PRINTM(conteiner.ukfCaSph->correctStruct.P);        
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