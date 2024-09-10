// #pragma once
// #include "ukf.h"
// #include "imm.h"

// template<class M,
//          template <typename> class StateModel,
//          template <typename> class MeasureModel,
//          template <typename> class ControlFunc>

// struct InitUKFStateModelCVMeasureModelSph
// {
//     M X0;
//     M P0;
//     M  processNoise;
//     M  measurementNoise;
//     ParamSigmaPoints p;
//     ControlFunc<M> controlFunc;

//     std::unique_ptr<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>> make_estimator()
//     {
//         return std::make_unique<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>>(X0, processNoise, measurementNoise, p); // создание IMM
//     }

//     void InitializationEstimator(const Detection<M>& detection)
//     {
//         typedef Eigen::SparseMatrix<double> SpMat;
//         typedef Eigen::Triplet<double> T;

//         // M Hp(3,6);
//         // Hp << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//         //       0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//         //       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;

//         SpMat Hp(3,6);
//         std::vector<T> tripletList;
//         tripletList.reserve(3);

//         tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 2, 1.0));
//         tripletList.push_back(T(2, 4, 1.0));
//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());
        
//         X0 = Hp.transpose() * Utils<M>::sph2CartMeas(detection.measurement);
        
//         //-------------------------------------------------------------------------
//         double process_var = 0.00001;
//         double sko_range  = 100.0;
//         double sko_Az = 0.1/3.0;
//         double sko_Um = 0.1/3.0;
//         double sko_Vr = 5.0;
//         p.alpha = 1e-3;
//         p.beta = 2.0;
//         p.kappa = 3.0 - X0.rows();
//         //-------------------------------------------------------------------------
     
//         processNoise.resize(3,3);
//         processNoise <<  process_var,            0.0,          0.0,
//                                 0.0,        process_var,       0.0,
//                                 0.0,             0.0,      process_var;

        
//         if (detection.measurement.rows() == 3)
//         {
//             measurementNoise.resize(3,3);
//             measurementNoise << pow(sko_range,2),          0.0,                  0.0,
//                                         0.0,         pow(sko_Az,2),              0.0,
//                                         0.0,                0.0,            pow(sko_Um,2);
//         }
//         else
//         {
//             measurementNoise.resize(4,4);
//             measurementNoise << pow(sko_range,2),          0.0,              0.0,          0.0,
//                                         0.0,         pow(sko_Az,2),          0.0,          0.0,
//                                         0.0,               0.0,         pow(sko_Um,2),     0.0,
//                                         0.0,               0.0,              0.0,        pow(sko_Vr,2);
//         }


//     }
// };

// template<class M,
//          template <typename> class StateModel,
//          template <typename> class MeasureModel,
//          template <typename> class ControlFunc>

// struct InitUKFStateModelCTMeasureModelSph
// {
//     M X0;
//     M P0;
//     M  processNoise;
//     M  measurementNoise;
//     ParamSigmaPoints p;
//     ControlFunc<M> controlFunc;

//     std::unique_ptr<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>> make_estimator()
//     {
//         return std::make_unique<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>>(X0, processNoise, measurementNoise, p);
//     }

//     void InitializationEstimator(const Detection<M>& detection)
//     {
//         typedef Eigen::SparseMatrix<double> SpMat;
//         typedef Eigen::Triplet<double> T;

//         // M Hp(3,7);
//         // Hp << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//         //       0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//         //       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;

//         SpMat Hp(3,7);
//         std::vector<T> tripletList;
//         tripletList.reserve(3);

//         tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 2, 1.0));
//         tripletList.push_back(T(2, 4, 1.0));
//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());
        
//         X0 = Hp.transpose() * Utils<M>::sph2CartMeas(detection.measurement);
        
//         //-------------------------------------------------------------------------
//         double process_var = 10.0;
//         double sko_range  = 100.0;
//         double sko_Az = 0.1/3.0;
//         double sko_Um = 0.1/3.0;
//         p.alpha = 1e-3;
//         p.beta = 2.0;
//         p.kappa = 3.0 - X0.rows();
//         //-------------------------------------------------------------------------
     
//         processNoise.resize(4,4);
//         processNoise <<  process_var,            0.0,          0.0,     0.0,
//                                 0.0,        process_var,       0.0,     0.0,
//                                 0.0,             0.0,          process_var,     0.0,
//                                 0.0,             0.0,          0.0,    1e-7;

        
//         measurementNoise.resize(3,3);
//         measurementNoise << pow(sko_range,2),          0.0,                  0.0,
//                                     0.0,         pow(sko_Az,2),              0.0,
//                                     0.0,                0.0,            pow(sko_Um,2);

//     }
// };


// template<class M,
//          template <typename> class StateModel,
//          template <typename> class MeasureModel,
//          template <typename> class ControlFunc>

// struct InitUKFStateModelCAMeasureModelSph
// {
//     M X0;
//     M P0;
//     M  processNoise;
//     M  measurementNoise;
//     ParamSigmaPoints p;
//     ControlFunc<M> controlFunc;

//     std::unique_ptr<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>> make_estimator()
//     {
//         return std::make_unique<UnscentedKalmanfilter<M, StateModel, MeasureModel, ControlFunc>>(X0, processNoise, measurementNoise, p);
//     }

//     void InitializationEstimator(const Detection<M>& detection)
//     {
//         typedef Eigen::SparseMatrix<double> SpMat;
//         typedef Eigen::Triplet<double> T;


//         SpMat Hp(3,9);
//         std::vector<T> tripletList;
//         tripletList.reserve(3);

//         tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 3, 1.0));
//         tripletList.push_back(T(2, 6, 1.0));
//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());
        
//         X0 = Hp.transpose() * Utils<M>::sph2CartMeas(detection.measurement);

        
//         //-------------------------------------------------------------------------
//         double process_var = 10.0;
//         double sko_range  = 100.0;
//         double sko_Az = 0.1/3.0;
//         double sko_Um = 0.1/3.0;
//         p.alpha = 1e-3;
//         p.beta = 2.0;
//         p.kappa = 3.0 - X0.rows();
//         //-------------------------------------------------------------------------
     
//         processNoise.resize(3,3);
//         processNoise <<  process_var,            0.0,          0.0,
//                                 0.0,        process_var,       0.0,
//                                 0.0,             0.0,      process_var;

        
//         measurementNoise.resize(3,3);
//         measurementNoise << pow(sko_range,2),          0.0,                  0.0,
//                                     0.0,         pow(sko_Az,2),              0.0,
//                                     0.0,                0.0,            pow(sko_Um,2);


                
//     }
// };


// template<class M,
//          template <typename> class ConteinerType>        
// struct InitUkfImmMeasureModelSph
// {
//     ConteinerType<M> conteiner; 
    
//     M  mu_i;
//     M  p_ij;
//     ParamSigmaPoints param;

//     std::unique_ptr<Imm<M, ConteinerType>> make_estimator()
//     {
//         return std::make_unique<Imm<M, ConteinerType>>(mu_i, p_ij, conteiner);
//     }

//     void InitializationEstimator(const Detection<M>& detection)
//     {
//         conteiner.initConteiner(detection.measurement);
//         mu_i.resize(1,3);     
//         mu_i << 1.0/3.0, 1.0/3.0, 1.0/3.0;

//         p_ij.resize(3,3);                   
//         p_ij << 0.97, 0.015, 0.015,
//                 0.015, 0.97, 0.015,
//                 0.015, 0.015, 0.97;
                                                                        
//     }

// };

// template<class M, class TypeEstimator, class TypeEstimatorInit>
// struct Track
// {
// private:
//     double timePoint;
// public:
//     std::unique_ptr<TypeEstimator> estimator;
//     Track(const Detection<M>& detection)
//     {
//         TypeEstimatorInit estimatorInit;
//         estimatorInit.InitializationEstimator(detection);
//         estimator = estimatorInit.make_estimator();
//         timePoint = detection.time;
//     }
//     M step(const Detection<M> &detection)
//     {
//         try
//         {   
//             double dt = detection.time - timePoint;
//             timePoint = detection.time;
//             return estimator->step(detection.measurement, dt);
//         }
//         catch (const std::runtime_error &e)
//         {
//             std::cerr << e.what() << '\n';
//             return M();
//         }
//     }
    
//     M step(double t)
//     {
//         try
//         {
//         double dt = t - timePoint;
//         timePoint = t;
//         return estimator->step(dt);
//         }
//                 catch (const std::exception &e)
//         {
//             std::cerr << e.what() << '\n';
//             return M();
//         }
//     }
// };
