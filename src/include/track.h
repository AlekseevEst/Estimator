#pragma once
#include "ukf.h"


template<class M,
         template <typename> class StateModel,
         template <typename> class MeasureModel,
         template <typename> class ControlFunc>

struct InitUKFStateModelCAMeasureModelSph
{
    M X0;
    M P0;
    M  processNoise;
    M  measurementNoise;
    ParamSigmaPoints p;
    ControlFunc<M> controlFunc;

    std::unique_ptr<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>> make_estimator()
    {
        return std::make_unique<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>>(X0, processNoise, measurementNoise, p);
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
     
        M Qp (3,3);
        Qp <<  process_var,            0.0,          0.0,
                                0.0,        process_var,       0.0,
                                0.0,             0.0,      process_var;


        M G = controlFunc(dt); 
        processNoise = G * Qp * G.transpose();

        
        measurementNoise.resize(3,3);
        measurementNoise << pow(sko_range,2),          0.0,                  0.0,
                                    0.0,         pow(sko_Az,2),              0.0,
                                    0.0,                0.0,            pow(sko_Um,2);

    }
};


template<class M, class TypeEstimator, class TypeEstimatorInit>
struct Track
{
private:
    std::unique_ptr<TypeEstimator> estimator;
    double timePoint;
public:
    Track(const Detection<M>& detection)
    {
        TypeEstimatorInit estimatorInit;
        estimatorInit.InitializationEstimator(detection);
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
        estimator->correctStruct.X = estimator->predict(dt);
        estimator->correctStruct.P = estimator->predictStruct.Pe;
        return estimator->correctStruct.X;
        }
                catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
            return M();
        }
    }
};
