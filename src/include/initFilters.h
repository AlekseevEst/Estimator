# pragma once
#include "ifilter.h"
#include "imm.h"
#include "ukf.h"
#include "utils.h"

struct InitUnscentedKalmanFilterCV
{
    template <class M,
              class... TypeArgs>
    void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCV, TypeArgs...> &filter,
                    const M &meas,
                    const M &measNoise)
    {
        M Hp(3, 6);
        Hp << 1, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 1, 0;

        filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);


        auto cartCov = Utils<M>::sph2cartcov(measNoise, meas);
        M posCov = cartCov.first;
        M velCov = cartCov.second;


        M Hv(3, 6);
        Hv << 0, 1, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1;

        filter.correctInfo.P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv;

        // filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());

        filter.Q << 0.00001, 0.0, 0.0,
                    0.0, 0.00001, 0.0,
                    0.0, 0.0, 0.00001;

        filter.R = measNoise;

        filter.paramsSigmaPoints.alpha = 1e-3;
        filter.paramsSigmaPoints.beta = 2.0;
        filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
    }
};

struct InitUnscentedKalmanFilterCT
{
    template <class M,
              class... TypeArgs>
    void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCT, TypeArgs...> &filter,
                    const M &meas,
                    const M &measNoise)

    {

        M Hp(3, 7);
        Hp << 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0;

        filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);

        auto cartCov = Utils<M>::sph2cartcov(measNoise, meas);
        M posCov = cartCov.first;
        M velCov = cartCov.second;

        M Hv(3, 7);
        Hv << 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0;

        M Hw(3, 7);
        Hw << 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1;

        M omega(3, 3);
        omega << 0, 0, 0,
                 0, 0, 0,
                 0, 0, pow(22.,2);

        filter.correctInfo.P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv + Hw.transpose() * omega * Hw;

        filter.Q.resize(4, 4); 
        filter.Q << 10.0, 0.0, 0.0, 0.0,
                     0.0, 10.0, 0.0, 0.0,
                     0.0, 0.0, 10.0, 0.0,
                     0.0, 0.0, 0.0, 1e-7;

        filter.R = measNoise;

        filter.paramsSigmaPoints.alpha = 1e-3;
        filter.paramsSigmaPoints.beta = 2.0;
        filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
    }
};

struct InitUnscentedKalmanFilterCA
{
    template <class M,
              class... TypeArgs>
    void operator()(UnscentedKalmanFilter<M, InitUnscentedKalmanFilterCA, TypeArgs...> &filter,
                    const M &meas,
                    const M &measNoise)
    {

        M Hp(3,9);
        Hp << 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0;

        filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);

        auto cartCov = Utils<M>::sph2cartcov(measNoise, meas);
        M posCov = cartCov.first;
        M velCov = cartCov.second;


        M Hv(3,9);
        Hv << 0, 1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0;

        M Ha(3,9);
        Ha << 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1;

        M AccelerationCov = M::Zero(3,3);
        AccelerationCov.diagonal() << pow(50,2), pow(50,2), pow(50,2);

        filter.correctInfo.P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv + Ha.transpose() * AccelerationCov * Ha;


        // filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());
        filter.Q.resize(3,3);
        filter.Q << 10.0, 0.0, 0.0,
                    0.0, 10.0, 0.0,
                    0.0, 0.0, 10.0;

        filter.R = measNoise;

        filter.paramsSigmaPoints.alpha = 1e-3;
        filter.paramsSigmaPoints.beta = 2.0;
        filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
    }
};


template <class M>
struct InitImmFilter1
{
    struct InitIMMFilterCV
    {
        template <class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterCV, TypeArgs...> &filter,
                        const M &meas,
                        const M &measNoise)
        {
           M Hp(3, 6);
        Hp << 1, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 1, 0;

        filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);


        auto cartCov = Utils<M>::sph2cartcov(measNoise, meas);
        M posCov = cartCov.first;
        M velCov = cartCov.second;


        M Hv(3, 6);
        Hv << 0, 1, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1;

        filter.correctInfo.P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv;

            // filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());

            filter.Q << 0.00001, 0.0, 0.0,
                        0.0, 0.00001, 0.0,
                        0.0, 0.0, 0.00001;

            filter.R = measNoise;

            filter.paramsSigmaPoints.alpha = 1e-3;
            filter.paramsSigmaPoints.beta = 2.0;
            filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
            // инициализируем фильтр
        }
    };

    using TypeFilterCV = UnscentedKalmanFilter<M, InitIMMFilterCV, FuncConstVel<M>, FuncMeasSph<M>, FuncControlMatrix_XvXYvYZvZ<M>>;

    struct InitIMMFilterCTxy
    {
        template <class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterCTxy, TypeArgs...> &filter,
                        const M &meas,
                        const M &measNoise)
        {
            M Hp(3, 7);
            Hp << 1, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0;

            filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);

            auto cartCov = Utils<M>::sph2cartcov(measNoise, meas);
            M posCov = cartCov.first;
            M velCov = cartCov.second;

            M Hv(3, 7);
            Hv << 0, 1, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 1, 0;

            M Hw(3, 7);
            Hw << 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1;

            M omega(3, 3);
            omega << 0, 0, 0,
                     0, 0, 0,
                     0, 0, pow(22,2);

            filter.correctInfo.P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv + Hw.transpose() * omega * Hw;

            // filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());
            filter.Q.resize(4, 4);
            filter.Q << 10.0, 0.0, 0.0, 0.0,
                        0.0, 10.0, 0.0, 0.0,
                        0.0, 0.0, 10.0, 0.0,
                        0.0, 0.0, 0.0, 1e-7;

            filter.R = measNoise;

            filter.paramsSigmaPoints.alpha = 1e-3;
            filter.paramsSigmaPoints.beta = 2.0;
            filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
        }
    };
    using TypeFilterCTxy = UnscentedKalmanFilter<M, InitIMMFilterCTxy, FuncConstTurnXY<M>, FuncMeasSph<M>, FuncControlMatrix_XvXYvYZvZW<M>>;

    struct InitIMMFilterCA
    {
        template <class... TypeArgs>
        void operator()(UnscentedKalmanFilter<M, InitIMMFilterCA, TypeArgs...> &filter,
                        const M &meas,
                        const M &measNoise)
        {
            M Hp(3, 9);
            Hp << 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0, 0;

            filter.correctInfo.X = Hp.transpose() * Utils<M>::sph2CartMeas(meas);

            auto cartCov = Utils<M>::sph2cartcov(measNoise, meas);
            M posCov = cartCov.first;
            M velCov = cartCov.second;

            M Hv(3, 9);
            Hv << 0, 1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0;

            M Ha(3, 9);
            Ha << 0, 0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 1;


            M AccelerationCov = M::Zero(3, 3);
            AccelerationCov.diagonal() << pow(50, 2), pow(50, 2), pow(50, 2);


            filter.correctInfo.P = Hp.transpose() * posCov * Hp + Hv.transpose() * velCov * Hv + Ha.transpose() * AccelerationCov * Ha;

            // filter.correctInfo.P = Utils<M>::do_cart_P0(Utils<M>::sph2cartcov(measNoise, meas), filter.correctInfo.X.rows());

            filter.Q << 10.0, 0.0, 0.0,
                        0.0, 10.0, 0.0,
                        0.0, 0.0, 10.0;

            filter.R = measNoise;

            filter.paramsSigmaPoints.alpha = 1e-3;
            filter.paramsSigmaPoints.beta = 2.0;
            filter.paramsSigmaPoints.kappa = 3.0 - filter.correctInfo.X.rows();
        }
    };
    using TypeFilterCA = UnscentedKalmanFilter<M, InitIMMFilterCA, FuncConstAcceleration<M>, FuncMeasSph<M>, FuncControlMatrix_XvXaXYvYaYZvZaZ<M>>;

    std::vector<std::shared_ptr<IFilter<M>>> filters; // <- если нужно проитерироваться по фильтрам

    std::shared_ptr<TypeFilterCV> filterCV; //<- для других моделей также
    std::shared_ptr<TypeFilterCTxy> filterCTxy;
    std::shared_ptr<TypeFilterCA> filterCA;
    /*по аналогии для других моделей*/

    InitImmFilter1() : filterCV{std::make_shared<TypeFilterCV>()},
                       filterCTxy{std::make_shared<TypeFilterCTxy>()},
                       filterCA{std::make_shared<TypeFilterCA>()}
    {
        filters.push_back(filterCV);
        filters.push_back(filterCTxy);
        filters.push_back(filterCA);
    }

    void operator()(IMM<M, InitImmFilter1> & filterIMM,
                    const M &meas,
                    const M &measNoise)
                    
    {
        filterCV->Initialization(meas, measNoise);
        filterCTxy->Initialization(meas,measNoise);
        filterCA->Initialization(meas,measNoise);
        filterIMM.mu_i << 1./3., 1./3., 1./3.;
        filterIMM.p_ij << 0.97, 0.015, 0.015,
                          0.015, 0.97, 0.015,
                          0.015, 0.015, 0.97;

    }
};

