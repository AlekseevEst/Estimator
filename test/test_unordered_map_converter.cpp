
// #include <catch2/catch.hpp>
// #include <iostream>
// #include <typeinfo>
// #include <typeindex>
// #include <unordered_map>
// #include <tuple>
// #include "keyMap.h"
// #include "models.h"
// #include <Eigen/Sparse>
// using namespace Catch::Benchmark;

// TEST_CASE("unordered_map")
// {

//     using SpMat = Eigen::SparseMatrix<double>;
//     using T = Eigen::Triplet<double>;

//     SpMat HVelTurn(6, 7);
//     SpMat HVelAcc(6, 9);
//     SpMat HTurnAcc(7, 9);

//     std::vector<T> tripletList;
//     tripletList.reserve(6);
//     tripletList.push_back(T(0, 0, 1)); // использовать POS_X
//     tripletList.push_back(T(1, 1, 1));
//     tripletList.push_back(T(2, 2, 1));
//     tripletList.push_back(T(3, 3, 1));
//     tripletList.push_back(T(4, 4, 1));
//     tripletList.push_back(T(5, 5, 1));
//     HVelTurn.setFromTriplets(tripletList.begin(), tripletList.end());
//     tripletList.clear();

//     tripletList.push_back(T(0, 0, 1));
//     tripletList.push_back(T(1, 1, 1));
//     tripletList.push_back(T(2, 3, 1));
//     tripletList.push_back(T(3, 4, 1));
//     tripletList.push_back(T(4, 6, 1));
//     tripletList.push_back(T(5, 7, 1));
//     HVelAcc.setFromTriplets(tripletList.begin(), tripletList.end());
//     HTurnAcc.setFromTriplets(tripletList.begin(),tripletList.end());
//     tripletList.clear();

//     FuncConstVel<Eigen::MatrixXd> modelCv;
//     FuncConstAcceleration<Eigen::MatrixXd> modelCa;
//     FuncConstTurn<Eigen::MatrixXd> modelCt;

//     std::unordered_map<KeyMap, std::function<Eigen::MatrixXd(Eigen::MatrixXd)>, Hash_fn> m;

//     m[{typeid(modelCv), typeid(modelCt)}] = [HVelTurn](Eigen::MatrixXd matStateOrCov)
//     {
//         Eigen::MatrixXd res;

//         if (matStateOrCov.cols() == 1)
//         {
//             res = HVelTurn.transpose() * matStateOrCov;
//             return res;

//         }
//         res = HVelTurn.transpose() * matStateOrCov * HVelTurn;
//         return res;

//     };

//     m[{typeid(modelCv), typeid(modelCa)}] = [HVelAcc](Eigen::MatrixXd matStateOrCov)
//     {
//         Eigen::MatrixXd res;

//         if (matStateOrCov.cols() == 1)
//         {
//             res = HVelAcc.transpose() * matStateOrCov;
//             return res;
//         }
//         res = HVelAcc.transpose() * matStateOrCov * HVelAcc;
//         return res;
//     };

//         m[{typeid(modelCt), typeid(modelCa)}] = [HTurnAcc](Eigen::MatrixXd matStateOrCov)
//     {
//         Eigen::MatrixXd res;

//         if (matStateOrCov.cols() == 1)
//         {
//             res = HTurnAcc.transpose() * matStateOrCov;
//             return res;
//         }
//         res = HTurnAcc.transpose() * matStateOrCov * HTurnAcc;
//         return res;
//     };

       
//         m[{typeid(modelCa), typeid(modelCv)}] = [HVelAcc](Eigen::MatrixXd matStateOrCov)
//     {
//         Eigen::MatrixXd res;

//         if (matStateOrCov.cols() == 1)
//         {
//             res = HVelAcc * matStateOrCov;
//             return res;
//         }
//         res = HVelAcc * matStateOrCov * HVelAcc.transpose();
//         return res;
//     };





//     if (m.count({typeid(modelCv), typeid(modelCt)}))
//     {
//         Eigen::MatrixXd exampleState(6, 1);
//         Eigen::MatrixXd exampleCov(6, 6);
//         exampleState << 100, 100, 100, 100, 100, 100;

//         exampleCov << 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100;

//         // std::cout << m[{typeid(modelCv), typeid(modelCt)}](exampleState) << std::endl;
//         // std::cout << m[{typeid(modelCv), typeid(modelCt)}](exampleCov) << std::endl;
//     }

//     if (m.count({typeid(modelCv), typeid(modelCa)}))
//     {
//         Eigen::MatrixXd exampleState(6, 1);
//         Eigen::MatrixXd exampleCov(6, 6);
//         exampleState << 100, 100, 100, 100, 100, 100;

//         exampleCov << 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100;

//         exampleState << 100, 100, 100, 100, 100, 100;
//         // std::cout << m[{typeid(modelCv), typeid(modelCa)}](exampleState) << std::endl;
//         // std::cout << m[{typeid(modelCv), typeid(modelCa)}](exampleCov) << std::endl;
//     }


//     if (m.count({typeid(modelCt), typeid(modelCa)}))
//     {
//         Eigen::MatrixXd exampleState(7, 1);
//         Eigen::MatrixXd exampleCov(7, 7);
//         exampleState << 100, 100, 100, 100, 100, 100, 100;

//         exampleCov << 100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100, 100;

//         // std::cout << m[{typeid(modelCt), typeid(modelCa)}](exampleState) << std::endl;
//         // std::cout << m[{typeid(modelCt), typeid(modelCa)}](exampleCov) << std::endl;
//     }

//     if (m.count({typeid(modelCa), typeid(modelCv)}))
//     {
//         Eigen::MatrixXd exampleState(9, 1);
//         Eigen::MatrixXd exampleCov(9, 9);
//         exampleState << 100, 100, 100, 100, 100, 100, 100, 100, 100;

//         exampleCov << 100, 100, 100, 100, 100, 100, 100, 100 ,100,
//                     100, 100, 100, 100, 100, 100, 100,  100 ,100,
//                     100, 100, 100, 100, 100, 100, 100, 100 ,100,
//                     100, 100, 100, 100, 100, 100, 100, 100 ,100,
//                     100, 100, 100, 100, 100, 100, 100, 100 ,100,
//                     100, 100, 100, 100, 100, 100, 100, 100 ,100,
//                     100, 100, 100, 100, 100, 100, 100, 100 ,100,
//                     100 ,100, 100 ,100, 100 ,100, 100 ,100, 100,
//                     100 ,100, 100 ,100, 100 ,100, 100 ,100, 100;

//         std::cout << m[{typeid(modelCa), typeid(modelCv)}](exampleState) << std::endl;
//         std::cout << m[{typeid(modelCa), typeid(modelCv)}](exampleCov) << std::endl;
//     }

//     BENCHMARK("unordered_map"){

//     };
// }
