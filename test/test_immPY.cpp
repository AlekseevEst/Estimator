// #include <iostream>
// #include <catch2/catch.hpp>
// #include "imm.h"
// #include "ukf.h"
// // #include "ifilter.h"
// // #include "initFilters.h"
// using namespace Catch::Benchmark;


// template <class M>
// struct InitImmFilterTest
// {
//     struct InitIMMFilterCV
//     {
//         template <class... TypeArgs>
//         void operator()(UnscentedKalmanFilter<M, InitIMMFilterCV, TypeArgs...> &filter,
//                         const M &meas,
//                         const M &measNoise)
//         {

//             filter.correctInfo.X << 18., 1. ,18., 1., 18.,1.;
            
//             filter.correctInfo.P << 1., 0, 0, 0, 0, 0,
//                                     0, 1. ,0 ,0 ,0 ,0,
//                                     0, 0, 1. ,0 , 0 ,0,
//                                     0, 0, 0 ,1. , 0 ,0,
//                                     0, 0, 0 ,0 , 1. ,0,
//                                     0, 0, 0 ,0 , 0 ,1.;

//             filter.Q.resize(3,3);
//             filter.Q << 1.0, 0.0, 0.0,
//                         0.0, 1.0, 0.0,
//                         0.0, 0.0, 1.0;

//             filter.R = measNoise;

//             filter.paramsSigmaPoints.alpha = 1e-3;
//             filter.paramsSigmaPoints.beta = 2.0;
//             filter.paramsSigmaPoints.kappa = 0.0;
//             // инициализируем фильтр
//         }
//     };

//     using TypeFilterCV = UnscentedKalmanFilter<M, InitIMMFilterCV, FuncConstVel<M>, FuncMeasSph<M>, FuncControlMatrix_XvXYvYZvZ<M>>;
 
//     struct InitIMMFilterCV2
//     {
//         template <class... TypeArgs>
//         void operator()(UnscentedKalmanFilter<M, InitIMMFilterCV2, TypeArgs...> &filter,
//                         const M &meas,
//                         const M &measNoise)
//         {

//             filter.correctInfo.X << 20., 2. ,20., 2., 20.,2.;
            
//             filter.correctInfo.P << 2., 0, 0, 0, 0, 0,
//                                     0, 2. ,0 ,0 ,0 ,0,
//                                     0, 0, 2. ,0 , 0 ,0,
//                                     0, 0, 0 ,2. , 0 ,0,
//                                     0, 0, 0 ,0 , 2. ,0,
//                                     0, 0, 0 ,0 , 0 ,2.;

//             filter.Q.resize(3,3);
//             filter.Q << 2.0, 0.0, 0.0,
//                         0.0, 2.0, 0.0,
//                         0.0, 0.0, 2.0;

//             filter.R = measNoise;

//             filter.paramsSigmaPoints.alpha = 1e-3;
//             filter.paramsSigmaPoints.beta = 2.0;
//             filter.paramsSigmaPoints.kappa = 0.0;
//             // инициализируем фильтр
//         }
//     };
//     using TypeFilterCV2 = UnscentedKalmanFilter<M, InitIMMFilterCV2, FuncConstVel<M>, FuncMeasSph<M>, FuncControlMatrix_XvXYvYZvZ<M>>;


//         struct InitIMMFilterCV3
//     {
//         template <class... TypeArgs>
//         void operator()(UnscentedKalmanFilter<M, InitIMMFilterCV3, TypeArgs...> &filter,
//                         const M &meas,
//                         const M &measNoise)
//         {

//             filter.correctInfo.X << 22., 3. ,22., 3., 22.,3.;
            
//             filter.correctInfo.P << 3., 0, 0, 0, 0, 0,
//                                     0, 3. ,0 ,0 ,0 ,0,
//                                     0, 0, 3. ,0 , 0 ,0,
//                                     0, 0, 0 ,3. , 0 ,0,
//                                     0, 0, 0 ,0 , 3. ,0,
//                                     0, 0, 0 ,0 , 0 ,3.;

//             filter.Q.resize(3,3);
//             filter.Q << 3.0, 0.0, 0.0,
//                         0.0, 3.0, 0.0,
//                         0.0, 0.0, 3.0;

//             filter.R = measNoise;

//             filter.paramsSigmaPoints.alpha = 1e-3;
//             filter.paramsSigmaPoints.beta = 2.0;
//             filter.paramsSigmaPoints.kappa = 0.0;
//             // инициализируем фильтр
//         }
//     };
//     using TypeFilterCV3 = UnscentedKalmanFilter<M, InitIMMFilterCV3, FuncConstVel<M>, FuncMeasSph<M>, FuncControlMatrix_XvXYvYZvZ<M>>;

//     std::vector<std::shared_ptr<IFilter<M>>> filters; // <- если нужно проитерироваться по фильтрам

//     std::shared_ptr<TypeFilterCV> filterCV;
//     std::shared_ptr<TypeFilterCV2> filterCV2;
//     std::shared_ptr<TypeFilterCV3> filterCV3; //<- для других моделей также

  
//     /*по аналогии для других моделей*/

//     InitImmFilterTest() : filterCV{std::make_shared<TypeFilterCV>()},
//                           filterCV2{std::make_shared<TypeFilterCV2>()},
//                           filterCV3{std::make_shared<TypeFilterCV3>()}


//     {
//         filters.push_back(filterCV);
//         filters.push_back(filterCV2);
//         filters.push_back(filterCV3);

//     }

//     void operator()(IMM<M, InitImmFilterTest> & filterIMM,
//                     const M &meas,
//                     const M &measNoise)
                    
//     {
//         filterCV->Initialization(meas, measNoise);
//         filterCV2->Initialization(meas,measNoise);
//         filterCV3->Initialization(meas,measNoise);


//         filterIMM.mu_i << 0.5,0.25,0.25;

//         filterIMM.p_ij << 0.97,0.015,0.015,
//                           0.015,0.97,0.015,
//                           0.015 ,0.015,0.97;
           
//     }
// };

// TEST_CASE("imm_step")
// {
//     Eigen::MatrixXd meas(3,1);
//     Eigen::MatrixXd measNoise(3,3);
//     measNoise << 1.,0, 0,
//                  0, 1.,0,
//                  0, 0, 1.;
 


//     IMM<Eigen::MatrixXd,InitImmFilterTest<Eigen::MatrixXd>> imm;
//     imm.Initialization(meas, measNoise);
//     double dt = 0.2;
                  
//     Eigen::MatrixXd expectedXpred(6,1);
//     expectedXpred << 19.82635843,  1.73925383, 19.82635843,  1.73925383, 19.82635843,  1.73925383;

//     Eigen::MatrixXd expectedPpred(6,6);
//     expectedPpred <<5.12779649, 1.86315645, 3.3182725,  1.50830568, 3.3182725,  1.50830568,
//                     1.86315645, 2.49484732, 1.50830568, 0.68559349, 1.50830568, 0.68559349,
//                     3.3182725,  1.50830568, 5.12779648, 1.86315645, 3.3182725,  1.50830568,
//                     1.50830568, 0.68559349, 1.86315645, 2.49484732, 1.50830568, 0.68559349,
//                     3.3182725,  1.50830568, 3.3182725,  1.50830568, 5.12779648, 1.86315645,
//                     1.50830568, 0.68559349, 1.50830568, 0.68559349, 1.86315645, 2.49484732;


//     Eigen::MatrixXd expectedXcor(6,1);
//     expectedXcor << 18.33072555,  1.73380428, 18.33072554,  1.73380428, 24.88462677,  3.01925687;
//     Eigen::MatrixXd expectedPcor(6,6);
//     expectedPcor << 0.44442207, 0.11517105, 0.19771506, 0.06678267, 0.16121659, 0.05963365,
//                     0.11517105, 2.44654534, 0.06678267, 0.07144681, 0.13733088, 0.08532779,
//                     0.19771506, 0.06678267, 0.44442207, 0.11517105, 0.16121658, 0.05963365,
//                     0.06678267, 0.07144681, 0.11517105, 2.44654534, 0.13733087, 0.08532779,
//                     0.16121659, 0.13733088, 0.16121658, 0.13733087, 0.7075415,  0.24459036,
//                     0.05963365, 0.08532779, 0.05963365, 0.08532779, 0.24459036, 2.48727504;
                 

//     Eigen::MatrixXd expectedUpdateModeProb(1,3);
//     expectedUpdateModeProb << 0.00579003, 0.58558497, 0.408625;

//     auto pred = imm.predict(dt);
//     CHECK((pred.first - expectedXpred).norm() == Approx(0.0).margin(1e-5));
//     bool condition = (pred.first - expectedXpred).norm() == Approx(0.0).margin(1e-5);
//     if(!condition){
//         PRINTM(pred.first.format(customFormat));
//         PRINTM(expectedXpred.format(customFormat));
//     }

//     CHECK((pred.second - expectedPpred).norm() == Approx(0.0).margin(1e-5));
//     condition = (pred.second - expectedPpred).norm() == Approx(0.0).margin(1e-5);
//     if(!condition){
//         PRINTM(pred.second.format(customFormat));
//         PRINTM(expectedPpred.format(customFormat));
//     }


//         // PRINTM(imm.initializator.filterCV->getPredictInfo().Xe.format(customFormat));
//         // PRINTM(imm.initializator.filterCV2->getPredictInfo().Xe.format(customFormat));
//         // PRINTM(imm.initializator.filterCV3->getPredictInfo().Xe.format(customFormat));
//         // PRINTM(imm.mu_ij.format(customFormat));



//     Eigen::MatrixXd Z (3,1);
//     Z << 35.355, 45., 45.; 
    
//     auto cor = imm.correct(Z);

//     CHECK((cor.first - expectedXcor).norm() == Approx(0.0).margin(1e-1));
//     condition = (cor.first - expectedXcor).norm() == Approx(0.0).margin(1e-1);
//     if(!condition){
//         PRINTM(cor.first.format(customFormat));
//         PRINTM(expectedXcor.format(customFormat));
//     }

//     CHECK((cor.second - expectedPcor).norm() == Approx(0.0).margin(1e-1));
//     condition = (cor.second - expectedPcor).norm() == Approx(0.0).margin(1e-1);
//     if(!condition){
//         PRINTM(cor.second.format(customFormat));
//         PRINTM(expectedPcor.format(customFormat));
//     }
//     CHECK((imm.mu_i - expectedUpdateModeProb).norm() == Approx(0.0).margin(1e-4));
//     condition = (imm.mu_i - expectedUpdateModeProb).norm() == Approx(0.0).margin(1e-4);
//     if(!condition){
//         PRINTM(imm.mu_i.format(customFormat));
//         PRINTM(expectedUpdateModeProb.format(customFormat));
//     }


//     BENCHMARK("imm_step"){
//     // auto pred = imm.predict(dt);
//     // auto cor = imm.correct(Z);
//     };

// }
    

   