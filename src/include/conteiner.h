#pragma once
#include "utils.h"
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
        ukfCvSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstVel, FuncMeasSph, FuncControlMatrix_XvXYvYZvZ>>(X0, procNoise, measNoise, paramSigmaPoints);


        procNoise.resize (4, 4);
        procNoise <<    10.0, 0.0, 0.0, 0.0,
                         0.0, 10.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 0.0, 0.0, 1e-7;

        paramSigmaPoints.kappa = -4.0;
        ukfCtSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstTurnXY, FuncMeasSph, FuncControlMatrix_XvXYvYZvZW>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCtXy)}](X0), procNoise, measNoise, paramSigmaPoints);

        process_var = 10.0;
        procNoise.resize(3, 3);
        procNoise <<    process_var, 0.0, 0.0,
                        0.0, process_var, 0.0,
                        0.0, 0.0, process_var;

        paramSigmaPoints.kappa = -6.0;
        ukfCaSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstAcceleration, FuncMeasSph, FuncControlMatrix_XvXaXYvYaYZvZaZ>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](X0), procNoise, measNoise, paramSigmaPoints);

        filters.push_back(ukfCvSph);
        filters.push_back(ukfCtSph);
        filters.push_back(ukfCaSph);
    }
};

template<class M>
struct ConteinerCVCACTxz
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
        ukfCvSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstVel, FuncMeasSph, FuncControlMatrix_XvXYvYZvZ>>(X0, procNoise, measNoise, paramSigmaPoints);


        procNoise.resize (4, 4);
        procNoise <<    10.0, 0.0, 0.0, 0.0,
                         0.0, 10.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 0.0, 0.0, 1e-7;

        paramSigmaPoints.kappa = -4.0;
        ukfCtSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstTurnXZ, FuncMeasSph, FuncControlMatrix_XvXYvYZvZW>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCtXy)}](X0), procNoise, measNoise, paramSigmaPoints);

        process_var = 10.0;
        procNoise.resize(3, 3);
        procNoise <<    process_var, 0.0, 0.0,
                        0.0, process_var, 0.0,
                        0.0, 0.0, process_var;

        paramSigmaPoints.kappa = -6.0;
        ukfCaSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstAcceleration, FuncMeasSph, FuncControlMatrix_XvXaXYvYaYZvZaZ>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](X0), procNoise, measNoise, paramSigmaPoints);

        filters.push_back(ukfCvSph);
        filters.push_back(ukfCtSph);
        filters.push_back(ukfCaSph);
    }

};