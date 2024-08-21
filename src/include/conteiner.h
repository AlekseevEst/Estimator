#pragma once
#include "utils.h"
template<class M>
struct ConteinerCVCACTxy
{
    Converter<M> converter;
    std::shared_ptr <IFilter<M>> ukfCvSph;
    std::shared_ptr <IFilter<M>> ukfCtSph;
    std::shared_ptr <IFilter<M>> ukfCaSph;

    ConteinerCVCACTxy(const M& X, const M& procNoise, const M& measNoise, ParamSigmaPoints paramSigmaPoints)
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
        ukfCtSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstTurnXY, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCt)}](X), ProcNoiseCT, measNoise, paramSigmaPoints);
        paramSigmaPoints.kappa = 3 - 9; //параметр ансцентного преобразования задал вручную для CA. 
        ukfCaSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](X), procNoise, measNoise, paramSigmaPoints);

    }
    // initConteinerCVCACT()
    // {
        // нужно
    //     // инициализиция фильтров своим состоянием.  Чтобы не передавать в imm (P) и тд, а только p_ij
    // }
};

template<class M>
struct ConteinerCVCACTxz
{
    Converter<M> converter;
    std::shared_ptr <IFilter<M>> ukfCvSph;
    std::shared_ptr <IFilter<M>> ukfCtSph;
    std::shared_ptr <IFilter<M>> ukfCaSph;

    ConteinerCVCACTxz(const M& X, const M& procNoise, const M& measNoise, ParamSigmaPoints paramSigmaPoints)
    {
        M ProcNoiseCT(4, 4);
        ProcNoiseCT <<   10.0, 0.0, 0.0, 0.0,
                         0.0, 10.0, 0.0, 0.0,
                         0.0, 0.0, 1.0, 0.0,
                         0.0, 0.0, 0.0, 1e-7;

        M ProcNoiseCV(3, 3);
        ProcNoiseCV <<   1.0, 0.0, 0.0,
                         0.0, 1.0, 0.0,
                         0.0, 0.0, 1.0;


        // здесь выделяем память. Пока что все вместе: и выделение и инициализация фильтров
        ukfCvSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>>(X, ProcNoiseCV, measNoise, paramSigmaPoints);
        paramSigmaPoints.kappa = 3 - 7; //параметр ансцентного преобразования задал вручную для CT. 
        ukfCtSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstTurnXZ, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCt)}](X), ProcNoiseCT, measNoise, paramSigmaPoints);
        paramSigmaPoints.kappa = 3 - 9; //параметр ансцентного преобразования задал вручную для CA. 
        ukfCaSph = std::make_shared<UnscentedKalmanfilter<M, FuncConstAcceleration, FuncMeasSphCA, FuncControlMatrix_XvXaXYvYaYZvZaZ>>(converter.m[{typeid(converter.modelCv), typeid(converter.modelCa)}](X), procNoise, measNoise, paramSigmaPoints);

    }
    // initConteinerCVCACT()
    // {
        // нужно
    //     // инициализиция фильтров своим состоянием.  Чтобы не передавать в imm (P) и тд, а только p_ij
    // }
};