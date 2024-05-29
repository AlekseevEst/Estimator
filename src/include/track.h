#pragma once
#include "ukf.h"


template<class M, class StateModel, class MeasureModel>
struct InitUKFStateModelCAMeasureModelSph
{
    M X0;
    M measurementNoise;
    M proccesNoise;
    ParamSigmaPoints p;

    std::unique_ptr<UnscentedKalmanfilter<M,StateModel,MeasureModel>> make_estimator()
    {
        return std::make_unique<UnscentedKalmanfilter<M,StateModel,MeasureModel>>(X0, proccesNoise, measurementNoise, p);
    }

    void InitializationEstimator(const Detection<M>& detection)
    {
        typedef Eigen::SparseMatrix<double> SpMat;
        typedef Eigen::Triplet<double> T;

        // M Hp(3,9);
        // Hp << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        //       0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        //       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;


        SpMat Hp(3,9);
        std::vector<T> tripletList;
        tripletList.reserve(3);

        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 3, 1.0));
        tripletList.push_back(T(2, 6, 1.0));
        Hp.setFromTriplets(tripletList.begin(), tripletList.end());
        
        X0 = Hp.transpose() * detection.point;
        
        //-------------------------------------------------------------------------
        double process_var = 10.0;
        double sko_range  = 100.0;
        double sko_Az = 0.1/3.0;
        double sko_Um = 0.1/3.0;
        double dt = 6.0;
        p.alpha = 1e-3;
        p.beta = 2.0;
        p.kappa = 3.0 - X0.rows();
        //-------------------------------------------------------------------------
     
        proccesNoise.resize (3,3);
        proccesNoise <<  process_var,            0.0,          0.0,
                                0.0,        process_var,       0.0,
                                0.0,             0.0,      process_var;
        
        measurementNoise.resize(3,3);
        measurementNoise << pow(sko_range,2),          0.0,                  0.0,
                                    0.0,         pow(sko_Az,2),              0.0,
                                    0.0,                0.0,            pow(sko_Um,2);

         
    }
};


template<class M, class EstimatorType, class EstimatorInit>
struct Track
{
private:
    std::unique_ptr<EstimatorType> estimator;
    double timePoint;
public:
    Track(const Detection<M>& detection)
    {
        EstimatorInit estimatorInit(detection);
        estimator = estimatorInit.make_estimator();
        timePoint = detection.timePoint;
    }
    M step(const Detection<M> &detection)
    {
        try
        {
            double dt = detection.timePoint - timePoint;
            timePoint = detection.timePoint;
            M Xe = estimator->predict(dt);
            M X = estimator->correct(detection.point);
            return X;
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
        estimator.correctStruct.X = estimator.predict(dt);
        estimator.correctStruct.P = estimator.predictStruct.Pe;
        return estimator.correctStruct.X;
        }
                catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return M();
        }
    }
};
