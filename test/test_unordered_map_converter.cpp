
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

//         SpMat HVelTurn(6, 7);
//         SpMat HVelAcc(6, 9);
//         SpMat HTurnAcc(7, 9);
//         SpMat HBalAcc(7,9);
//         SpMat HVelBal(6,7);
//         SpMat HTurnBal(7,7);

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

//             tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 1, 1.0));
//         tripletList.push_back(T(2, 3, 1.0));
//         tripletList.push_back(T(3, 4, 1.0));
//         tripletList.push_back(T(5, 6, 1.0));
//         tripletList.push_back(T(6, 7, 1.0));
//         HBalAcc.setFromTriplets(tripletList.begin(), tripletList.end());
//         tripletList.clear();

//         tripletList.push_back(T(0, 0, 1.0));
//         tripletList.push_back(T(1, 1, 1.0));
//         tripletList.push_back(T(2, 2, 1.0));
//         tripletList.push_back(T(3, 3, 1.0));
//         tripletList.push_back(T(4, 5, 1.0));
//         tripletList.push_back(T(5, 6, 1.0));
//         HVelBal.setFromTriplets(tripletList.begin(), tripletList.end());
//         HTurnBal.setFromTriplets(tripletList.begin(), tripletList.end());
//         tripletList.clear();

//     FuncConstVel<Eigen::MatrixXd> modelCv;
//     FuncConstAcceleration<Eigen::MatrixXd> modelCa;
//     FuncConstTurnXY<Eigen::MatrixXd> modelCt;
//     FuncBalreentry<Eigen::MatrixXd> modelBal;

//     std::unordered_map<namespaceKeyMap::KeyMap, std::function<Eigen::MatrixXd(Eigen::MatrixXd)>, namespaceKeyMap::Hash_fn> m;

//     // m[{typeid(modelCv), typeid(modelCt)}] = [HVelTurn](Eigen::MatrixXd matStateOrCov)
//     // {
//     //     Eigen::MatrixXd res;

//     //     if (matStateOrCov.cols() == 1)
//     //     {
//     //         res = HVelTurn.transpose() * matStateOrCov;
//     //         return res;

//     //     }
//     //     res = HVelTurn.transpose() * matStateOrCov * HVelTurn;
//     //     return res;

//     // };

//     // m[{typeid(modelCv), typeid(modelCa)}] = [HVelAcc](Eigen::MatrixXd matStateOrCov)
//     // {
//     //     Eigen::MatrixXd res;

//     //     if (matStateOrCov.cols() == 1)
//     //     {
//     //         res = HVelAcc.transpose() * matStateOrCov;
//     //         return res;
//     //     }
//     //     res = HVelAcc.transpose() * matStateOrCov * HVelAcc;
//     //     return res;
//     // };

//     //     m[{typeid(modelCt), typeid(modelCa)}] = [HTurnAcc](Eigen::MatrixXd matStateOrCov)
//     // {
//     //     Eigen::MatrixXd res;

//     //     if (matStateOrCov.cols() == 1)
//     //     {
//     //         res = HTurnAcc.transpose() * matStateOrCov;
//     //         return res;
//     //     }
//     //     res = HTurnAcc.transpose() * matStateOrCov * HTurnAcc;
//     //     return res;
//     // };

       
//     //     m[{typeid(modelCa), typeid(modelCv)}] = [HVelAcc](Eigen::MatrixXd matStateOrCov)
//     // {
//     //     Eigen::MatrixXd res;

//     //     if (matStateOrCov.cols() == 1)
//     //     {
//     //         res = HVelAcc * matStateOrCov;
//     //         return res;
//     //     }
//     //     res = HVelAcc * matStateOrCov * HVelAcc.transpose();
//     //     return res;
//     // };


//         m[{typeid(modelBal), typeid(modelCa)}] = [HBalAcc](Eigen::MatrixXd matStateOrCov)
//     {
//         Eigen::MatrixXd res;

//         if (matStateOrCov.cols() == 1)
//         {
//             res = HBalAcc.transpose() * matStateOrCov;
//             return res;
//         }
//         res = HBalAcc.transpose() * matStateOrCov * HBalAcc;
//         return res;
//     };

       
//         m[{typeid(modelCa), typeid(modelCv)}] = [HBalAcc](Eigen::MatrixXd matStateOrCov)
//     {
//         Eigen::MatrixXd res;

//         if (matStateOrCov.cols() == 1)
//         {
//             res = HBalAcc * matStateOrCov;
//             return res;
//         }
//         res = HBalAcc * matStateOrCov * HBalAcc.transpose();
//         return res;
//     };
//         m[{typeid(modelCv), typeid(modelBal)}] = [HVelBal](Eigen::MatrixXd matStateOrCov)
//     {
//         Eigen::MatrixXd res;

//             if (matStateOrCov.cols() == 1)
//             {
//                 res = HVelBal.transpose() * matStateOrCov;
//                 return res;
//             }
//             res = HVelBal.transpose() * matStateOrCov * HVelBal;
//             return res;
//     };



//     // if (m.count({typeid(modelCv), typeid(modelCt)}))
//     // {
//     //     Eigen::MatrixXd exampleState(6, 1);
//     //     Eigen::MatrixXd exampleCov(6, 6);
//     //     exampleState << 100, 100, 100, 100, 100, 100;

//     //     exampleCov << 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100;

//     //     // std::cout << m[{typeid(modelCv), typeid(modelCt)}](exampleState) << std::endl;
//     //     // std::cout << m[{typeid(modelCv), typeid(modelCt)}](exampleCov) << std::endl;
//     // }

//     // if (m.count({typeid(modelCv), typeid(modelCa)}))
//     // {
//     //     Eigen::MatrixXd exampleState(6, 1);
//     //     Eigen::MatrixXd exampleCov(6, 6);
//     //     exampleState << 100, 100, 100, 100, 100, 100;

//     //     exampleCov << 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100;

//     //     exampleState << 100, 100, 100, 100, 100, 100;
//     //     // std::cout << m[{typeid(modelCv), typeid(modelCa)}](exampleState) << std::endl;
//     //     // std::cout << m[{typeid(modelCv), typeid(modelCa)}](exampleCov) << std::endl;
//     // }


//     // if (m.count({typeid(modelCt), typeid(modelCa)}))
//     // {
//     //     Eigen::MatrixXd exampleState(7, 1);
//     //     Eigen::MatrixXd exampleCov(7, 7);
//     //     exampleState << 100, 100, 100, 100, 100, 100, 100;

//     //     exampleCov << 100, 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100, 100,
//     //                 100, 100, 100, 100, 100, 100, 100;

//     //     // std::cout << m[{typeid(modelCt), typeid(modelCa)}](exampleState) << std::endl;
//     //     // std::cout << m[{typeid(modelCt), typeid(modelCa)}](exampleCov) << std::endl;
//     // }

//     // if (m.count({typeid(modelCa), typeid(modelCv)}))
//     // {
//     //     Eigen::MatrixXd exampleState(9, 1);
//     //     Eigen::MatrixXd exampleCov(9, 9);
//     //     exampleState << 100, 100, 100, 100, 100, 100, 100, 100, 100;

//     //     exampleCov << 100, 100, 100, 100, 100, 100, 100, 100 ,100,
//     //                 100, 100, 100, 100, 100, 100, 100,  100 ,100,
//     //                 100, 100, 100, 100, 100, 100, 100, 100 ,100,
//     //                 100, 100, 100, 100, 100, 100, 100, 100 ,100,
//     //                 100, 100, 100, 100, 100, 100, 100, 100 ,100,
//     //                 100, 100, 100, 100, 100, 100, 100, 100 ,100,
//     //                 100, 100, 100, 100, 100, 100, 100, 100 ,100,
//     //                 100 ,100, 100 ,100, 100 ,100, 100 ,100, 100,
//     //                 100 ,100, 100 ,100, 100 ,100, 100 ,100, 100;

//     //     std::cout << m[{typeid(modelCa), typeid(modelCv)}](exampleState) << std::endl;
//     //     std::cout << m[{typeid(modelCa), typeid(modelCv)}](exampleCov) << std::endl;
//     // }


//         if (m.count({typeid(modelCa), typeid(modelBal)}))
//     {   
//         PRINT(1111111111);
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

//         std::cout << m[{typeid(modelCa), typeid(modelBal)}](exampleState) << std::endl;
//         std::cout << m[{typeid(modelCa), typeid(modelBal)}](exampleCov) << std::endl;
//     }



//     if (m.count({typeid(modelCv), typeid(modelBal)}))
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
//         std::cout << m[{typeid(modelCv), typeid(modelBal)}](exampleState) << std::endl;
//         std::cout << m[{typeid(modelCv), typeid(modelBal)}](exampleCov) << std::endl;
//     }


//         if (m.count({typeid(modelBal), typeid(modelCa)}))
//     {
//         Eigen::MatrixXd exampleState(7, 1);
//         Eigen::MatrixXd exampleCov(7, 7);

//         exampleCov << 100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100, 100,
//                     100, 100, 100, 100, 100, 100,100,
//                     100, 100, 100, 100, 100, 100,100,
//                     100, 100, 100, 100, 100, 100,100,
//                     100, 100, 100, 100, 100, 100,100,
//                     100, 100, 100, 100, 100, 100,100;

//         exampleState << 100, 100, 100, 100, 100, 100,100;
//         std::cout << m[{typeid(modelBal), typeid(modelCa)}](exampleState) << std::endl;
//         std::cout << m[{typeid(modelBal), typeid(modelCa)}](exampleCov) << std::endl;
//     }
//     BENCHMARK("unordered_map"){

//     };
// }
