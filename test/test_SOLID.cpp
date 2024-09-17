// #include <catch2/catch.hpp>
// #include <iostream>
// #include <typeinfo>
// #include <typeindex>

// #include <cstddef>

// #include "Eigen/Dense"
// #include "Eigen/Sparse"
// #include <cmath>
// #include <vector>
// #include <iostream>
// #include "utils.h"
// #include "models.h"
// #include "converter.h"

// using namespace Catch::Benchmark;

// // template <class M>
// // struct Detection {
// //     double time;
// //     M Measurement;
// //     M MeasurementNoise;
// //     struct DetectionParams {
// //         /*
// //             ....
// //         */
// //     };
// //     DetectionParams params;
// // };

// template <class M,
//           class TypeEstimator>
// struct Track
// {
//     Track(): estimator{std::make_unique<TypeEstimator>()}
//     {
//         // трассы лучше получить заранее, положить их в контейнер и отуда забирать
//         // в конструкторе выделяем память под всё что нужно, тут нужно создать класс типа TypeEstimator
//     }

//     void Initialization(const Detection<M> &detection)
//     {
//         // вызвать, обрати внимание дальше тип Detection не идёт

//         estimator->Initialization(detection.measurement,
//                                  detection.measurementNoise);
//         timePoint = detection.time;
//         // PRINTM(estimator->correctInfo.X);
//         // PRINTM(estimator->correctInfo.P);
//     }

//     M step(const Detection<M> &detection)
//     {
//         try
//         {   

//             // std::cout<<detection.time;
//             double dt = detection.time - timePoint;
//             timePoint = detection.time;
//             auto Pred = estimator->predict(dt);
//             // PRINTM(detection.Measurement);
//             return estimator->correct(detection.measurement).first;
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
//             double dt = t - timePoint;
//             timePoint = t;
//             return estimator->predict(dt).first;
//         }
//         catch (const std::exception &e)
//         {
//             std::cerr << e.what() << '\n';
//             return M();
//         }
//     }

// private:
//     double timePoint;
//     std::unique_ptr<TypeEstimator> estimator;
// };

// template <class M>
// struct UnscentedKalmanFilterMath
// {

//     void compute_weights(ParamSigmaPoints &paramSigmaPoints, double &lamda, double &c, const int &dim_x, std::vector<double> &Wc, std::vector<double> &Wm)
//     {
//         lamda = pow(paramSigmaPoints.alpha, 2) * (dim_x + paramSigmaPoints.kappa) - dim_x;
//         c = 0.5 / (dim_x + lamda);
//         Wc.assign(2 * dim_x + 1, c);
//         Wm.assign(2 * dim_x + 1, c);
//         Wc[0] = lamda / (dim_x + lamda) + (1 - pow(paramSigmaPoints.alpha, 2) + paramSigmaPoints.beta);
//         Wm[0] = lamda / (dim_x + lamda);
//     }

//     void compute_sigma_points(const M &X, const M &P, ParamSigmaPoints &paramSigmaPoints, const double &lamda, M &Xu, M &U)
//     {
//         U = sqrt(lamda + X.rows()) * Utils<M>::sqrtMatSpectral(P);

//         Xu.col(0) = X;

//         for (size_t i = 0; i < X.rows(); i++)
//         {
//             Xu.col(i + 1) = X + U.col(i);
//         }
//         for (size_t i = 0; i < X.rows(); i++)
//         {
//             Xu.col(i + X.rows() + 1) = X - U.col(i);
//         }
//     }

//     void doExtrapolatedStateVector(const M &Xue, M &Xe, const std::vector<double>& Wm)
//     {
//         //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------
//         Xe.setZero();
        
        
//         for (int i = 0; i < Xue.cols(); i++)
//         {
//             Xe += Wm[i] * Xue.col(i);
//         }
//     }

//     void doCovMatExtrapolatedStateVector(const M &Xue, const M &Xe, const M &G, const M &Q, M &Pe, const std::vector<double>& Wc)
//     {
//         //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ
//         Pe.setZero();
//         for (int i = 0; i < Xue.cols(); i++)
//         {
          
//             Pe += Wc[i] * ((Xue.col(i) - Xe) * (Xue.col(i) - Xe).transpose());
//         }

//         Pe += G * Q * G.transpose();
   
//     }

//     void doExtrapolatedMeasVector(const M &Zue, M &Ze, const std::vector<double>& Wm)
//     {
//         //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.  ------------------
//         Ze.setZero();
//         for (int i = 0; i < Zue.cols(); i++)
//         {
//             Ze += Wm[i] * Zue.col(i);
//         }
//     }

//     void doCovMatExtrapolatedMeasVector(const M &Zue, const M &Ze, M &Pzz, const std::vector<double>& Wc)
//     {
//         //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА ИЗМЕРЕНИИ СФЕР.
//         Pzz.setZero();
//         for (int i = 0; i < Zue.cols(); i++)
//         {

//             Pzz += Wc[i] * ((Zue.col(i) - Ze) * (Zue.col(i) - Ze).transpose());
//         }
//     }

//     void doCovMatInnovation(const M &Pzz, const M &R_sph_deg, M &Se)
//     {

//         Se = Pzz + R_sph_deg; // innovation covariance

//     }

//     //----------------------------------------------------------------------

//     void calcGainFilter(const M &Xue, const M &Xe, const M &Zue, const M &Ze, const M &Se, M &Pxz, M &gainKalman, const std::vector<double>& Wc)
//     {
//         Pxz.setZero();

//         for (int i = 0; i < Zue.cols(); i++)
//         {

//             Pxz += Wc[i] * (Xue.col(i) - Xe) * (Zue.col(i) - Ze).transpose();
//         }
      
//         gainKalman = Pxz * Se.inverse();
        
//     }

//     void correctState(const M &Xe, const M &Z, const M &Ze, const M &K, M &X)

//     {
//         X = Xe + K * (Z - Ze);
     
//     }

//     void correctCov(const M &Pe, const M &K, const M &Se, M &P)
//     {
//         P = Pe - K * Se * K.transpose();
//         if (!Utils<M>::CheckingConditionsMat(P)) // проверка на симметричность, положительно определённость и не вырожденность
//             throw std::runtime_error("СheckingСonditionsMat ERROR");
//     }

// private:

// };


// template <class M>
// struct IFilter
// {
//     virtual std::pair<M, M> predict(double dt) = 0;
//     virtual std::pair<M, M> correct(const M &Z) = 0;


//     virtual double likelihood(/*...*/) = 0;
//     virtual double distance(const M &Z) = 0;

//     virtual std::type_index getModelType() const = 0;

//     virtual void setCorrectInfo(const M& X, const M& P) = 0;
//     virtual Predict<M> getPredictInfo() = 0;
//     virtual Correct<M> getCorrectInfo() = 0;

//     virtual ~IFilter() {}
// };

// template <class M,
//           class TypeInitialization,
//           class TypeStateModelFunc,
//           class TypeMeasFunc,
//           class TypeControlFunc>
// struct UnscentedKalmanFilter
//     : public IFilter<M>
// {
//     std::pair<M, M> predict(double dt) override final
//     {
//         UKfilterMath.compute_weights(paramsSigmaPoints, lamda, c, correctInfo.X.rows(), Wc, Wm);
//         UKfilterMath.compute_sigma_points(correctInfo.X, correctInfo.P, paramsSigmaPoints, lamda, sigmaVectors, U);
//         extrapolatedStateSigmaVectors = stateFunc(sigmaVectors, dt);
//         UKfilterMath.doExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictInfo.Xe, Wm);
//         G = controlFunc(dt);
//         UKfilterMath.doCovMatExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictInfo.Xe, G, Q, predictInfo.Pe, Wc);

//         correctInfo.X = predictInfo.Xe;  // Записываю в качестве скореектированных предсказанные значения. под вопросом! Нужно в track.step(dt)
//         correctInfo.P = predictInfo.Pe;

//         return std::make_pair(predictInfo.Xe, predictInfo.Pe);
//     }

//     std::pair<M, M> correct(const M &Z) override final
//     {
//         extrapolatedMeasSigmaVectors = measFunc(extrapolatedStateSigmaVectors, Z);
//         UKfilterMath.doExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictInfo.Ze, Wm);
//         UKfilterMath.doCovMatExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictInfo.Ze, predictInfo.Pzz, Wc);
//         UKfilterMath.doCovMatInnovation(predictInfo.Pzz, R, predictInfo.Se);
//         UKfilterMath.calcGainFilter(extrapolatedStateSigmaVectors, predictInfo.Xe, extrapolatedMeasSigmaVectors, predictInfo.Ze, predictInfo.Se, predictInfo.Pxz, predictInfo.K, Wc);
//         UKfilterMath.correctState(predictInfo.Xe, Z, predictInfo.Ze, predictInfo.K, correctInfo.X);
//         UKfilterMath.correctCov(predictInfo.Pe, predictInfo.K, predictInfo.Se, correctInfo.P);

//         return std::make_pair(correctInfo.X, correctInfo.P);
//     }

//     double likelihood(/*...*/) override final
//     {

//         // Что должно быть здесь? функция правдоподобия высчитывается, вроде бы только в IMM алгоритме.
//     }



//     double distance(const M &Z) override final
//     {

//         v = Z - predictInfo.Ze;
//         double distance = (v.transpose() * predictInfo.Se.inverse() * v)(0,0);
//         return distance;
//     }

//     std::type_index getModelType() const override
//     {
//         return std::type_index(typeid(stateFunc));
//     }

//     Correct<M> getCorrectInfo() override final {

//         return correctInfo;
//     }

//     Predict<M> getPredictInfo() override final {

//         return predictInfo;
//     }
//     void setCorrectInfo(const M& X, const M& P) override final
//     {
//         correctInfo.X = X;
//         correctInfo.P = P;
//     }

//     UnscentedKalmanFilter()
//     {
//         predictInfo.Xe.resize(stateFunc.getSize(), 1);
//         predictInfo.Pe.resize(stateFunc.getSize(), stateFunc.getSize());
//         predictInfo.Ze.resize(measFunc.getSize(), 1);
//         predictInfo.Pzz.resize(measFunc.getSize(), measFunc.getSize());
//         predictInfo.Pxz.resize(stateFunc.getSize(), measFunc.getSize());
//         predictInfo.Se.resize(measFunc.getSize(), measFunc.getSize());
//         predictInfo.K.resize(stateFunc.getSize(), measFunc.getSize());
//         correctInfo.X.resize(stateFunc.getSize(), 1);
//         correctInfo.P.resize(stateFunc.getSize(), stateFunc.getSize());
//         Q.resize(measFunc.getSize(), measFunc.getSize());
//         R.resize(measFunc.getSize(), measFunc.getSize());
//         G.resize(controlFunc.getSize().first, controlFunc.getSize().second);
//         v.resize(measFunc.getSize(), 1);
//         sigmaVectors.resize(stateFunc.getSize(), 2 * stateFunc.getSize() + 1);
//         Wc.reserve(sigmaVectors.cols());
//         Wm.reserve(sigmaVectors.cols());
//         U.resize(correctInfo.P.rows(), correctInfo.P.cols());
//         extrapolatedStateSigmaVectors.resize(sigmaVectors.rows(), sigmaVectors.cols());
//         extrapolatedMeasSigmaVectors.resize(measFunc.getSize(), extrapolatedStateSigmaVectors.cols());
//     }

//     template <class... TypeArgs>
//     void Initialization(TypeArgs... args)
//     {
//         initializator(*this, args...);
//     }

//     Predict<M> predictInfo;
//     Correct<M> correctInfo;
//     ParamSigmaPoints paramsSigmaPoints;
//     TypeStateModelFunc stateFunc;
//     TypeMeasFunc measFunc;
//     M Q;
//     M R;

// private:
//     TypeInitialization initializator;

//     TypeControlFunc controlFunc;
//     UnscentedKalmanFilterMath<M> UKfilterMath;
//     M extrapolatedStateSigmaVectors;
//     M extrapolatedMeasSigmaVectors;
//     M G;
//     M v;
//     M sigmaVectors;
//     M U;
//     double lamda;
//     double c;
//     std::vector<double> Wc, Wm;
// };

// struct InitUnscentedKalmanFilterCV
// {
//     template <class M,
//               class... TypeArgs>
//     void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCV, TypeArgs...> &filter,
//                     const M &meas,
//                     const M &measNoise)
//     {
//         using SpMat = Eigen::SparseMatrix<double>;
//         using T = Eigen::Triplet<double>;

//         SpMat Hp(3, 6);
//         std::vector<T> tripletList;
//         tripletList.reserve(3);

//         tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 2, 1.0));
//         tripletList.push_back(T(2, 4, 1.0));
//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//         filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);

//         filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());

//         filter.Q << 0.00001, 0.0, 0.0,
//                     0.0, 0.00001, 0.0,
//                     0.0, 0.0, 0.00001;

//         filter.R = measNoise;

//         filter.paramsSigmaPoints.alpha = 1e-3;
//         filter.paramsSigmaPoints.beta = 2.0;
//         filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
//     }
// };

// struct InitUnscentedKalmanFilterCT
// {
//     template <class M,
//               class... TypeArgs>
//     void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCT, TypeArgs...> &filter,
//                     const M &meas,
//                     const M &measNoise)
    

//         {
//             using SpMat = Eigen::SparseMatrix<double>;
//             using T = Eigen::Triplet<double>;

//             SpMat Hp(3, 7);
//             std::vector<T> tripletList;
//             tripletList.reserve(3);

//             tripletList.push_back(T(0, 0, 1.0));
//             tripletList.push_back(T(1, 2, 1.0));
//             tripletList.push_back(T(2, 4, 1.0));
//             Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//             filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);
//             filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());
//             filter.Q.resize(4,4); // !!!!!!!!!!!!!!!!!!!!!!!!!!!
//             filter.Q << 10.0, 0.0, 0.0, 0.0,
//                          0.0, 10.0, 0.0, 0.0,
//                          0.0, 0.0, 10.0, 0.0,
//                          0.0, 0.0, 0.0, 1e-7;

//             filter.R = measNoise;

//             filter.paramsSigmaPoints.alpha = 1e-3;
//             filter.paramsSigmaPoints.beta = 2.0;
//             filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
//         }
    
// };

// struct InitUnscentedKalmanFilterCA
// {
//     template <class M,
//               class... TypeArgs>
//     void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCA, TypeArgs...> &filter,
//                     const M &meas,
//                     const M &measNoise)
//     {

//         using SpMat = Eigen::SparseMatrix<double>;
//         using T = Eigen::Triplet<double>;

//         SpMat Hp(3, 9);
//         std::vector<T> tripletList;
//         tripletList.reserve(3);

//         tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 3, 1.0));
//         tripletList.push_back(T(2, 6, 1.0));
//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//         filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);
//         filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());

//         filter.Q << 10.0, 0.0, 0.0,
//                     0.0, 10.0, 0.0,
//                     0.0, 0.0, 10.0;

//         filter.R = measNoise;

//         filter.paramsSigmaPoints.alpha = 1e-3;
//         filter.paramsSigmaPoints.beta = 2.0;
//         filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
//     }
// };

// // UnscentedKalmanFilter<Eigen::MatrixXd,
// //                       InitUnscentedKalmanFilterCV,
// //                       FuncConstVel<Eigen::MatrixXd>,
// //                       FuncMeasSph<Eigen::MatrixXd>,
// //                       FuncControlMatrix_XvXYvYZvZ<Eigen::MatrixXd>> filterInitCV;
    

// // UnscentedKalmanFilter<Eigen::MatrixXd,
// //                       InitUnscentedKalmanFilterCT,
// //                       FuncConstTurnXY<Eigen::MatrixXd>,
// //                       FuncMeasSph<Eigen::MatrixXd>,
// //                       FuncControlMatrix_XvXYvYZvZW<Eigen::MatrixXd>> filterInitCTxy;
                      
// // UnscentedKalmanFilter<Eigen::MatrixXd,
// //                       InitUnscentedKalmanFilterCA,
// //                       FuncConstAcceleration<Eigen::MatrixXd>,
// //                       FuncMeasSph<Eigen::MatrixXd>,
// //                       FuncControlMatrix_XvXaXYvYaYZvZaZ<Eigen::MatrixXd>> filterInitCA;




// template <class M> 
// struct ImmMath{

//     void computeMixingProbability(const M &p_ij, const M &mu_i, M &mu_ij, M &cj)
//     {

//         for (long int j = 0; j < mu_i.cols(); j++)
//         {
//             double c = 0.;
//             for (long int i = 0; i < p_ij.cols(); i++)
//             {
//                 c +=(p_ij(i, j) * mu_i(0, i));
//             }
//             cj(0, j) = c;
//         }

//         for (long int j = 0; j < p_ij.cols(); j++)
//         {
//             for (long int i = 0; i < p_ij.rows(); i++)
//             {
//                 mu_ij(i, j) = p_ij(i, j) * mu_i(0, i) / cj(0, j);
//             }
//         }
//     }

//     void MixingStateAndCovariance(const M &mu_ij, std::vector<M> &stateMixed, std::vector<M> &covarianceMixed, const std::vector<std::shared_ptr<IFilter<M>>> &filters, Converter<M> &converter)
//     {
//         for (size_t j = 0; j < filters.size(); ++j)
//         {
//             stateMixed[j].setZero();
//             covarianceMixed[j].setZero();
//         }

//         for (size_t j = 0; j < filters.size(); ++j)
//         {
            
//             for (size_t i = 0; i < filters.size(); ++i)

//             {
//                 M convertedState = converter.m[{filters[i]->getModelType(), filters[j]->getModelType()}](filters[i]->getCorrectInfo().X);
//                 stateMixed[j] += mu_ij(i, j) * convertedState;
//             }
//             for (size_t i = 0; i < filters.size(); ++i)
//             {
//                 M convertedState = converter.m[{filters[i]->getModelType(),filters[j]->getModelType()}](filters[i]->getCorrectInfo().X);
//                 M dX = convertedState - stateMixed[j];
//                 M convertedCovariance = converter.m[{filters[i]->getModelType(), filters[j]->getModelType()}](filters[i]->getCorrectInfo().P);
//                 covarianceMixed[j] += mu_ij(i, j) * (convertedCovariance + dX * dX.transpose());
//             }
//         }
//     }

//     double likelihoodFunction(const M &Z, const M &Ze, const M &Se)
//     {
//         M v = Z - Ze;
//         long double power = -0.5 * (v.transpose() * Se.inverse() * v)(0, 0);
//         long double probability = std::pow((1 / (2 * M_PI)), Z.rows() / 2.0) / std::sqrt(Se.determinant()) * std::exp(power);

//         return probability;
//     }

//     void updateModeProbability(const M &Z, const M &cj, M& mu_i, std::vector<std::shared_ptr<IFilter<M>>>& filters)
//     {
//         long double c = 0.0;
//         for (size_t i = 0; i < filters.size(); ++i)
//         {
//             mu_i(0, i) = likelihoodFunction(Z, filters[i]->getPredictInfo().Ze, filters[i]->getPredictInfo().Se) * cj(0, i);
//             c += mu_i(0, i);
//         }

//         for (size_t i = 0; i < filters.size(); ++i)
//         {
//             mu_i(0, i) /= c;
//         }

//     }

//         std::pair<M,M> combinationModelsCondition(M& mu_i, M& X, M& P, std::vector<std::shared_ptr<IFilter<M>>>& filters, Converter<M>& converter) 
//     {
//         X.setZero();
//         P.setZero();
//         for (size_t i = 0; i < filters.size(); ++i)
//         {
//             M convertedState = converter.m[{filters[i]->getModelType(), typeid(converter.modelCv)}](filters[i]->getCorrectInfo().X);
//             X += mu_i(0, i) * convertedState;
//         }

//         for (size_t i = 0; i < filters.size(); ++i)
//         {
//             M dx = converter.m[{filters[i]->getModelType(), typeid(converter.modelCv)}](filters[i]->getCorrectInfo().X) - X;
//             M convertedCovariance = converter.m[{filters[i]->getModelType(), typeid(converter.modelCv)}](filters[i]->getCorrectInfo().P);
//             P += mu_i(0, i) * (convertedCovariance + dx * dx.transpose());
//         }

//         return std::make_pair(X,P);
//     }


// private:

// };

// template <class M,
//           class TypeInitialization> 
// struct IMM
//     : public IFilter<M>
// {
//     std::pair<M, M> predict(double dt) override final
//     {
//         math.computeMixingProbability(p_ij, mu_i, mu_ij, cj);

//         math.MixingStateAndCovariance(mu_ij, stateMixed, covarianceMixed, initializator.filters, converter);

//         for (size_t i = 0; i < initializator.filters.size(); ++i)
//             {
//                 initializator.filters[i]->setCorrectInfo(stateMixed[i],covarianceMixed[i]);
//                 // initializator.filters[i]->correctInfo.X = stateMixed[i]; //SET X
//                 // initializator.filters[i]->correctInfo.P = covarianceMixed[i]; // SET P
//                 initializator.filters[i]->predict(dt);
//                 mu_i(0, i) = cj(0, i);
//             }

//         return math.combinationModelsCondition(mu_i, predictInfo.Xe, predictInfo.Pe, initializator.filters, converter);
//     }
      
//     std::pair<M, M> correct(const M &Z) override final
//     {
//         for (size_t i = 0; i < initializator.filters.size(); ++i)
//         {
//             initializator.filters[i]->correct(Z);
//         }

//         math.updateModeProbability(Z, cj, mu_i, initializator.filters);
//         return math.combinationModelsCondition(mu_i, correctInfo.X, correctInfo.P, initializator.filters, converter);
//     }
    

//     double likelihood(/*...*/) override final
//     { //????
//     }

//      double distance(const M &Z) override final
//     {
//         double totalDistance;
//         for (size_t i = 0; i < initializator.filters.size(); ++i)
//         {
//         M v = Z - initializator.filters[i]->getPredictInfo().Ze;
//         double mahalonobisDistance = (v.transpose() * initializator.filters[i]->getPredictInfo().Se.inverse() * v)(0,0);
//         totalDistance += mahalonobisDistance * mu_i(0,i);
//         }
//         return totalDistance;
//     }

//     std::type_index getModelType() const override
//     {
//         //????
//     }

//     Correct<M> getCorrectInfo() override final {

//         return correctInfo;
//     }

//     Predict<M> getPredictInfo() override final {

//         return predictInfo;
//     }
//     void setCorrectInfo(const M& X, const M& P) override final
//     {
//         correctInfo.X = X;
//         correctInfo.P = P;
//     }

//     IMM()
//     {
//         mu_i.resize(1, initializator.filters.size());
//         p_ij.resize(initializator.filters.size(),initializator.filters.size());
//         mu_ij.resize(p_ij.rows(),p_ij.cols());
//         cj.resize(1, mu_i.cols());
//         stateMixed.resize(initializator.filters.size());
//         covarianceMixed.resize(initializator.filters.size());

//         correctInfo.X.resize(initializator.filterCV->correctInfo.X.rows(), initializator.filterCV->correctInfo.X.cols());
//         correctInfo.P.resize(correctInfo.X.rows(), correctInfo.X.rows());
//         predictInfo.Xe.resize(initializator.filterCV->correctInfo.X.rows(), initializator.filterCV->correctInfo.X.cols());
//         predictInfo.Pe.resize(correctInfo.X.rows(), correctInfo.X.rows());


//         for (size_t i = 0; i < initializator.filters.size(); ++i)
//         {
//             stateMixed[i].resize(initializator.filters[i]->getCorrectInfo().X.rows(),initializator.filters[i]->getCorrectInfo().X.cols());
//             covarianceMixed[i].resize(stateMixed[i].rows(),stateMixed[i].rows());
//         }    

//     }

//     template <class... TypeArgs>
//     void Initialization(TypeArgs... args)
//     {
//         initializator(*this, args...);
//     }

//     M mu_i; // Вероятности режима i
//     M p_ij; // переходная вероятность режима из i в j
//     Predict<M> predictInfo;
//     Correct<M> correctInfo;

//     std::vector<M> stateMixed;
//     std::vector<M> covarianceMixed;
// private:
//     TypeInitialization initializator;
//     ImmMath<M> math;
//     Converter<M> converter;

//     M mu_ij; // смешенная вероятность
//     M cj;


// };

// template <class M>
// struct InitImmFilter1
// {
//     struct InitIMMFilterCV
//     {
//         template <class... TypeArgs>
//         void operator()(UnscentedKalmanFilter<M, InitIMMFilterCV, TypeArgs...> &filter,
//                         const M &meas,
//                         const M &measNoise)
//         {
//             using SpMat = Eigen::SparseMatrix<double>;
//             using T = Eigen::Triplet<double>;

//             SpMat Hp(3, 6);
//             std::vector<T> tripletList;
//             tripletList.reserve(3);

//             tripletList.push_back(T(0, 0, 1.0));
//             tripletList.push_back(T(1, 2, 1.0));
//             tripletList.push_back(T(2, 4, 1.0));
//             Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//             filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);

//             filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());

//             filter.Q << 0.00001, 0.0, 0.0,
//                         0.0, 0.00001, 0.0,
//                         0.0, 0.0, 0.00001;

//             filter.R = measNoise;

//             filter.paramsSigmaPoints.alpha = 1e-3;
//             filter.paramsSigmaPoints.beta = 2.0;
//             filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
//             // инициализируем фильтр
//         }
//     };

//     using TypeFilterCV = UnscentedKalmanFilter<M, InitIMMFilterCV, FuncConstVel<M>, FuncMeasSph<M>, FuncControlMatrix_XvXYvYZvZ<M>>;

//     struct InitIMMFilterCTxy
//     {
//         template <class... TypeArgs>
//         void operator()(UnscentedKalmanFilter<M, InitIMMFilterCTxy, TypeArgs...> &filter,
//                         const M &meas,
//                         const M &measNoise)
//         {
//             using SpMat = Eigen::SparseMatrix<double>;
//             using T = Eigen::Triplet<double>;

//             SpMat Hp(3, 7);
//             std::vector<T> tripletList;
//             tripletList.reserve(3);

//             tripletList.push_back(T(0, 0, 1.0));
//             tripletList.push_back(T(1, 2, 1.0));
//             tripletList.push_back(T(2, 4, 1.0));
//             Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//             filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);
//             filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());
//             filter.Q.resize(4,4);
//             filter.Q << 10.0, 0.0, 0.0, 0.0,
//                 0.0, 10.0, 0.0, 0.0,
//                 0.0, 0.0, 10.0, 0.0,
//                 0.0, 0.0, 0.0, 1e-7;

//             filter.R = measNoise;

//             filter.paramsSigmaPoints.alpha = 1e-3;
//             filter.paramsSigmaPoints.beta = 2.0;
//             filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
//         }
//     };
//     using TypeFilterCTxy = UnscentedKalmanFilter<M, InitIMMFilterCTxy, FuncConstTurnXY<M>, FuncMeasSph<M>, FuncControlMatrix_XvXYvYZvZW<M>>;

//     struct InitIMMFilterCA
//     {
//         template <class... TypeArgs>
//         void operator()(UnscentedKalmanFilter<M, InitIMMFilterCA, TypeArgs...> &filter,
//                         const M &meas,
//                         const M &measNoise)
//         {
//             using SpMat = Eigen::SparseMatrix<double>;
//             using T = Eigen::Triplet<double>;

//          SpMat Hp(3,9);
//         std::vector<T> tripletList;
//         tripletList.reserve(3);

//         tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 3, 1.0));
//         tripletList.push_back(T(2, 6, 1.0));
//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//         filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);
//         filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());

//         filter.Q <<    10.0,    0.0,    0.0,  
//                         0.0,   10.0,    0.0,
//                         0.0,    0.0,   10.0;
                      


//         filter.R = measNoise;

//         filter.paramsSigmaPoints.alpha = 1e-3;
//         filter.paramsSigmaPoints.beta = 2.0;
//         filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
//         }
//     };
//     using TypeFilterCA = UnscentedKalmanFilter<M, InitIMMFilterCA, FuncConstAcceleration<M>, FuncMeasSph<M>, FuncControlMatrix_XvXaXYvYaYZvZaZ<M>>;

//     std::vector<std::shared_ptr<IFilter<M>>> filters; // <- если нужно проитерироваться по фильтрам

//     std::shared_ptr<TypeFilterCV> filterCV; //<- для других моделей также
//     std::shared_ptr<TypeFilterCTxy> filterCTxy;
//     std::shared_ptr<TypeFilterCA> filterCA;
//     /*по аналогии для других моделей*/

//     InitImmFilter1() : filterCV{std::make_shared<TypeFilterCV>()},
//                        filterCTxy{std::make_shared<TypeFilterCTxy>()},
//                        filterCA{std::make_shared<TypeFilterCA>()}
//     {
//         filters.push_back(filterCV);
//         filters.push_back(filterCTxy);
//         filters.push_back(filterCA);
//     }

//     void operator()(IMM<M, InitImmFilter1> & filterIMM,
//                     const M &meas,
//                     const M &measNoise)
                    
//     {
//         filterCV->Initialization(meas, measNoise);
//         filterCTxy->Initialization(meas,measNoise);
//         filterCA->Initialization(meas,measNoise);


//         filterIMM.mu_i << 1./3., 1./3., 1./3.;
//         filterIMM.p_ij << 0.97, 0.015, 0.015,
//                           0.015, 0.97, 0.015,
//                           0.015, 0.015, 0.97;


//     }
// };

// TEST_CASE("SOLID")
// {

//     Eigen::MatrixXd Z(3, 1);
//     Eigen::MatrixXd Z1(3, 1);
//     Eigen::MatrixXd Z2(3, 1);
//     Eigen::MatrixXd R(3, 3);

//     Z << 130000.6547, 0.001, 0.001;
//     Z1 << 131200.1248, 0.001, 0.001;    
//     Z2 << 132400.1248, 0.001, 0.001;          
        
//     R << 10000.0,         0.0,                0.0,          
//              0.0,    pow((0.1/3),2),          0.0,         
//              0.0,         0.0,           pow((0.1/3),2);
            


//     Detection<Eigen::MatrixXd> d;
//     d.time = 6.0;
//     d.measurement = Z;
//     d.measurementNoise = R;

//     Detection<Eigen::MatrixXd> d1;
//     d1.time = 12.0;
//     d1.measurement = Z1;
//     d1.measurementNoise = R;

//     Detection<Eigen::MatrixXd> d2;
//     d2.time = 18.0;
//     d2.measurement = Z2;
//     d2.measurementNoise = R;

//     Track<Eigen::MatrixXd, IMM<Eigen::MatrixXd, InitImmFilter1<Eigen::MatrixXd>>> track;
   
//     track.Initialization(d);
//     Eigen::MatrixXd res = track.step(d1);

//     PRINTM(res);
//     res = track.step(d2);
//     PRINTM(res);
//     res = track.step(24.0);
//     PRINTM(res);
//     BENCHMARK("STEP"){
    
//     };

// }
