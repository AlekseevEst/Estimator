#include <catch2/catch.hpp>
#include <iostream>
#include <typeinfo>
#include <typeindex>

#include <cstddef>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "immMath.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "utils.h"
#include "models.h"
#include "sigma_points.h"
using namespace Catch::Benchmark;

template <class M>
struct UnscentedKalmanFilterMath
{

    UnscentedKalmanFilterMath()
    {
    }

    M doSigmaVectors(const M &X, const M &P, ParamSigmaPoints paramSigmaPoints)
    {
        //----------СОЗДАЕМ Xu СИГМА-ВЕКТОРОВ------------------

        sigmaPoints.compute_weights(paramSigmaPoints);
        return sigmaPoints.compute_sigma_points(X, P, paramSigmaPoints);
    }

    M doExtrapolatedStateVector(const M &Xue)
    {
        //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ----------
        Xe.setZero(Xue.rows(), 1);
        for (int i = 0; i < Xue.cols(); i++)
        {
            // PRINTM(sigmaPoints.Wm[i]);
            Xe = Xe + sigmaPoints.Wm[i] * Xue.col(i);
        }
        // PRINTM(Xe);
        return Xe;
    }

    M doCovMatExtrapolatedStateVector(const M &Xue, const M &Xe, const M &G, const M &Q)
    {
        //-----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СОСТОЯНИЯ
        Pe.setZero(Xue.rows(), Xue.rows());

        for (int i = 0; i < Xue.cols(); i++)
        {
            dX = Xue.col(i) - Xe;
            Pe = Pe + sigmaPoints.Wc[i] * (dX * dX.transpose());
        }

        Pe = Pe + G * Q * G.transpose();
        // PRINTM(Pe);

        return Pe;
    }

    M doExtrapolatedMeasVector(const M &Zue)
    {
        //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА СФЕР.  ------------------
        Ze.setZero(Zue.rows(), 1);
        for (int i = 0; i < Zue.cols(); i++)
        {
            Ze = Ze + sigmaPoints.Wm[i] * Zue.col(i);
        }
        return Ze;
    }

    M doCovMatExtrapolatedMeasVector(const M &Zue, const M &Ze)
    {
        //----------СТАТИСТИЧЕСКАЯ ОЦЕНКА МАТРИЦЫ КОВАРИАЦИИ ЭКСТРАПОЛИРОВАННОГО ВЕКТОРА ИЗМЕРЕНИИ СФЕР.
        Pzz.setZero(Zue.rows(), Zue.rows());
        for (int i = 0; i < Zue.cols(); i++)
        {
            v = Zue.col(i) - Ze; // невязка
            Pzz = Pzz + sigmaPoints.Wc[i] * (v * v.transpose());
        }
        return Pzz;
    }

    M doCovMatInnovation(const M &Pzz, const M &R_sph_deg)
    {
        Se = Pzz + R_sph_deg; // innovation covariance
        // PRINTM(Se);
        return Se;
    }

    //----------------------------------------------------------------------

    M calcGainFilter(const M &Xue, const M &Xe, const M &Zue, const M &Ze, const M &Se)
    {
        M Pxz = M::Zero(Xue.rows(), Zue.rows());

        for (int i = 0; i < Zue.cols(); i++)
        {
            // PRINTM(sigmaPoints.Wc[i]);
            dX = Xue.col(i) - Xe;
            // PRINTM(dX);
            v = Zue.col(i) - Ze;
            // PRINTM(v);
            Pxz = Pxz + sigmaPoints.Wc[i] * dX * v.transpose();
        }
        // PRINTM(Pxz);
        gainKalman = Pxz * Se.inverse();
        // PRINTM(gainKalman);
        return gainKalman;
    }

    M correctState(const M &Xe, const M &Z, const M &Ze, const M &K)

    {
        X = Xe + K * (Z - Ze);
        // PRINTM(X);
        return X;
    }

    M correctCov(const M &Pe, const M &K, const M &Se)
    {
        P = Pe - K * Se * K.transpose();
        if ((Utils<M>::CheckingConditionsMat(P))) // проверка на симметричность, положительно определённость и не вырожденность
            return P;
        else
            throw std::runtime_error("СheckingСonditionsMat ERROR");
    }

private:
    M Xe;
    M Pe;
    M dX;
    M Ze;
    M Pzz;
    M v;
    M Se;
    M gainKalman;

    M X;
    M P;
    ParamSigmaPoints paramSigmaPoints;
    SigmaPoints<M> sigmaPoints;
};

// template <class M>
// struct Detection {
//     double Time;
//     M Measurement;
//     M MeasurementNoise;
//     struct DetectionParams {
//         /*
//             ....
//         */
//     };
//     DetectionParams params;
// };


template <class M,
          class TypeEstimator>
struct Track
{
    Track()
    {
        // трассы лучше получить заранее, положить их в контейнер и отуда забирать
        // в конструкторе выделяем память под всё что нужно, тут нужно создать класс типа TypeEstimator

        // Здесь я создаю объект типа TypeEstimator
        //  std::make_unique<TypeEstimator<....>>(...);
    }

    void Initialization(const Detection<M> &detection)
    {
        // вызвать, обрати внимание дальше тип Detection не идёт

        // Не до конца понял фразу выше. А куда оно дальше может поити?
        estimator.Initialization(detection.Measurement,
                                 detection.MeasurementNoise);
    }

    M Step(const Detection<M> &detection)
    {
    }

    M Step(double t)
    {
    }

private:
    double timePoint;
    std::unique_ptr<TypeEstimator> estimator;
};

template <class M>
struct IFilter
{
    virtual std::pair<M, M> predict(double dt) = 0;
    virtual std::pair<M, M> correct(const M &Z) = 0;

    virtual double likelihood(/*...*/) = 0;
    virtual double distance(const M &Z) = 0;

    virtual std::type_index getModelType() const = 0;
    virtual ~IFilter() {}
};

template <class M,
          class TypeInitialization,
          class TypeStateModelFunc,
          class TypeMeasFunc,
          class TypeControlFunc> 
struct UnscentedKalmanFilter
    : public IFilter<M>
{
    std::pair<M, M> predict(double dt) override final
    {

        sigmaVectors = UKfilterMath.doSigmaVectors(correctInfo.X, correctInfo.P, paramsSigmaPoints);
        extrapolatedStateSigmaVectors = stateFunc(sigmaVectors, dt);
        predictInfo.Xe = UKfilterMath.doExtrapolatedStateVector(extrapolatedStateSigmaVectors);
        G = controlFunc(dt);
        predictInfo.Pe = UKfilterMath.doCovMatExtrapolatedStateVector(extrapolatedStateSigmaVectors, predictInfo.Xe, G, Q);

        return std::make_pair(predictInfo.Xe, predictInfo.Pe);
    }

    std::pair<M, M> correct(const M &Z) override final
    {
        extrapolatedMeasSigmaVectors = measFunc(extrapolatedStateSigmaVectors, Z);
        predictInfo.Ze = UKfilterMath.doExtrapolatedMeasVector(extrapolatedMeasSigmaVectors);
        predictInfo.Pzz = UKfilterMath.doCovMatExtrapolatedMeasVector(extrapolatedMeasSigmaVectors, predictInfo.Ze);
        predictInfo.Se = UKfilterMath.doCovMatInnovation(predictInfo.Pzz, R);
        predictInfo.K = UKfilterMath.calcGainFilter(extrapolatedStateSigmaVectors, predictInfo.Xe, extrapolatedMeasSigmaVectors, predictInfo.Ze, predictInfo.Se);
        correctInfo.X = UKfilterMath.correctState(predictInfo.Xe, Z, predictInfo.Ze, predictInfo.K);
        correctInfo.P = UKfilterMath.correctCov(predictInfo.Pe, predictInfo.K, predictInfo.Se);

        return std::make_pair(correctInfo.X, correctInfo.P);
    }

    double likelihood(/*...*/) override final
    {

        // Что должно быть здесь? функция правдоподобия высчитывается, вроде бы только в IMM алгоритме.
    }

    double distance(const M &Z) override final
    {

        // v = Z - predictInfo.Ze;
        // double distance = v.transpose() * predictInfo.Se.inverse() + predictInfo.Se.determinant();
        // return distance;
    }

    std::type_index getModelType() const override
    {
        return std::type_index(typeid(stateFunc));
    }

    UnscentedKalmanFilter()
    {
        predictInfo.Xe.resize(stateFunc.getSize(), measFunc.getSize());
        predictInfo.Pe.resize(stateFunc.getSize(), stateFunc.getSize());
        predictInfo.Ze.resize(measFunc.getSize(), measFunc.getSize());
        predictInfo.Pzz.resize(measFunc.getSize(), measFunc.getSize());
        predictInfo.Se.resize(measFunc.getSize(), measFunc.getSize());
        predictInfo.K.resize(stateFunc.getSize(), measFunc.getSize());

        correctInfo.X.resize(stateFunc.getSize(), measFunc.getSize());
        correctInfo.P.resize(stateFunc.getSize(), stateFunc.getSize());

        Q.resize(measFunc.getSize(), measFunc.getSize());
        R.resize(measFunc.getSize(), measFunc.getSize());
        G.resize(controlFunc.getSize().first, controlFunc.getSize().second);
        v.resize(measFunc.getSize(), measFunc.getSize());

        sigmaVectors.resize(stateFunc.getSize(), 2 * stateFunc.getSize() + 1);
        extrapolatedStateSigmaVectors.resize(sigmaVectors.rows(), sigmaVectors.cols());
        extrapolatedMeasSigmaVectors.resize(measFunc.getSize(), extrapolatedStateSigmaVectors.cols());
        
    }

    template <class... TypeArgs>
    void Initialization(TypeArgs... args)
    {
        initializator(*this, args...);
    }

    Predict<M> predictInfo;
    Correct<M> correctInfo;
    M Q;
    M R;
    ParamSigmaPoints paramsSigmaPoints;

private:
    TypeInitialization initializator;
    UnscentedKalmanFilterMath<M> UKfilterMath;
    M extrapolatedStateSigmaVectors;
    M extrapolatedMeasSigmaVectors;
    TypeStateModelFunc stateFunc;
    TypeMeasFunc measFunc;
    TypeControlFunc controlFunc;
    M G;
    M v;
    M sigmaVectors;
};

struct InitUnscentedKalmanFilterCV
{
    template <class M,
              class... TypeArgs>
    void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCV, TypeArgs...> &filter,
                    const M &meas,
                    const M &measNoise)
    {
        using SpMat = Eigen::SparseMatrix<double>;
        using T = Eigen::Triplet<double>;
        // typedef Eigen::SparseMatrix<double> SpMat;
        // typedef Eigen::Triplet<double> T;


        SpMat Hp(3,6);
        std::vector<T> tripletList;
        tripletList.reserve(3);

        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 2, 1.0));
        tripletList.push_back(T(2, 4, 1.0));
        Hp.setFromTriplets(tripletList.begin(), tripletList.end());
        
        filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);


        
        filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas),filter.correctInfo.X.rows());

        filter.Q <<     0.00001,         0.0,         0.0,
                           0.0,        0.00001,       0.0,
                           0.0,          0.0,       0.00001;

        filter.R = measNoise;

        filter.paramsSigmaPoints.alpha = 1e-3;
        filter.paramsSigmaPoints.beta = 2.0;
        filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();


    }
};

// struct InitUnscentedKalmanFilterCT {
//     template<class M,
//              class... TypeArgs>
//     void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCT, TypeArgs...>& filter,
//                     const M& meas,
//                     const M& measNoise) {
//         // инициализируем фильтр
//     }

    
// };

UnscentedKalmanFilter<Eigen::MatrixXd,
                      InitUnscentedKalmanFilterCV,
                      FuncConstVel<Eigen::MatrixXd>,
                      FuncMeasSph<Eigen::MatrixXd>,
                      FuncControlMatrix_XvXYvYZvZ<Eigen::MatrixXd>> filterInitCV;

// UnscentedKalmanFilter<Eigen::MatrixXd,
//                       InitUnscentedKalmanFilterCT,
//                       std::nullptr_t,
//                       std::nullptr_t,
//                       std::nullptr_t> filterInit2;


template <class M,
          class TypeInitialization>
struct IMM
        : public IFilter<M>
{
    std::pair<M, M> Predict(double dt) override final {

    }

    std::pair<M, M> Correct(const M &Z) override final {
    }

    double Likelihood(/*...*/) override final {
    }

    double Distance(/*...*/) override final {
    }

    std::type_index GetModelType() const override {
    }

    IMM() {
        //Выделяешь память тут
    }

    template <class ... TypeArgs>
    void Initialization(TypeArgs... args) {
        initializator(*this, args...);
    }

private:
    TypeInitialization initializator;
    //IMMmath<M> math;
    //.. и прочее

};

template <class M>
struct InitImmFilter1 {
    struct InitIMMFilterCV {
        template<class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterCV, TypeArgs...>& filter,
                        const M& meas,
                        const M& noise) {
                           
            // инициализируем фильтр
        }
    };

    using TypeFilterCV = UnscentedKalmanFilter<M, InitIMMFilterCV, /**/std::nullptr_t, std::nullptr_t,std::nullptr_t>;

    std::shared_ptr<TypeFilterCV> filterCV; //<- для других моделей также
    /*по аналогии для других моделей*/

    std::vector<std::shared_ptr<IFilter<M>>> filters; // <- если нужно проитерироваться по фильтрам

    InitImmFilter1() : filterCV{std::make_shared<TypeFilterCV>()} {
        filters.push_back(filterCV);
    }


    void operator()(IMM<M, InitImmFilter1>&,
                    const M& meas,
                    const M& noise) {
        filterCV->Initialization(meas, noise);

        // инициализируем фильтры
    }
};

TEST_CASE("SOLID") {

}

