#pragma once
#include <iostream>
#include "structs.h"

template <class M,
          class TypeEstimator>
struct Track
{
    Track(): estimator{std::make_unique<TypeEstimator>()}
    {
        // трассы лучше получить заранее, положить их в контейнер и отуда забирать
        // в конструкторе выделяем память под всё что нужно, тут нужно создать класс типа TypeEstimator
    }

    void Initialization(const Detection<M> &detection)
    {
        // вызвать, обрати внимание дальше тип Detection не идёт

        estimator->Initialization(detection.measurement,
                                 detection.measurementNoise);
        timePoint = detection.time;

    }

    M step(const Detection<M> &detection)
    {
        try
        {   

            double dt = detection.time - timePoint;
            timePoint = detection.time;
            auto Pred = estimator->predict(dt);
            return estimator->correct(detection.measurement).first;
        }
        catch (const std::runtime_error &e)
        {
            std::cerr << e.what() << '\n';
            return M();
        }
    }

    M step(double t)
    {
        try
        {
            double dt = t - timePoint;
            timePoint = t;
            return estimator->predict(dt).first;
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return M();
        }
    }

private:
    double timePoint;
    std::unique_ptr<TypeEstimator> estimator;
};