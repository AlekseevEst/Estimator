// #include <iostream>
// #include <Eigen/Dense>
// #include <catch2/catch.hpp>

// using namespace Catch::Benchmark;

// // template <class M>
// // Eigen::MatrixXd sqrtMat(const M& P)
// // {
// //     Eigen::LLT<M> lltofP(P);
// //     if (lltofP.info() != Eigen::Success)
// //     {
// //         throw std::runtime_error("Cholesky decomposition ERROR");
// //     }
// //     Eigen::MatrixXd L = lltofP.matrixL();
// //     return L;
// // }

// // TEST_CASE("c")
// // {
// //     Eigen::MatrixXd p(2, 2);
// //     p << 4,1,
// //         1,3;
// //     // p << 9999.94, 0, 17.134, 0, 17.1341, 0, 0,
// //     //      0, 40000, 0, 0, 0, 0, 0,
// //     //      17.134, 0, 182.906, 0, 0.0299046, 0, 0,
// //     //      0, 0, 0, 40000, 0, 0, 0,
// //     //      17.1341, 0, 0.0299046, 0, 182.907, 0, 0,
// //     //      0, 0, 0, 0, 0, 40000, 0,
// //     //      0, 0, 0, 0, 0, 0, 0;

// //     // Check eigenvalues
// //     Eigen::EigenSolver<Eigen::MatrixXd> es(p);
// //     std::cout << "Eigenvalues:\n" << es.eigenvalues().real() << std::endl;

// //     try {
// //         Eigen::MatrixXd res = sqrtMat(p);
// //         std::cout << res << std::endl;
// //         std::cout << res * res << std::endl;

// //     } catch (const std::exception& e) {
// //         std::cout << e.what() << std::endl;
// //     }

// //     BENCHMARK("bench"){

// //     };
// // }

// #include <iostream>
// #include <catch2/catch.hpp>
// #include "ukf.h"
// #include "converter.h"
// #include "ifilter.h"

// using namespace Catch::Benchmark;

// Eigen::MatrixXd sqrtMatSpectral(const Eigen::MatrixXd& P)
// {
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(P);
//     if (eigensolver.info() != Eigen::Success)
//     {
//         throw std::runtime_error("Eigen decomposition ERROR");
//     }

//     Eigen::MatrixXd D = eigensolver.eigenvalues().asDiagonal();
//     PRINTM(D);
//     Eigen::MatrixXd V = eigensolver.eigenvectors();

//     for (int i = 0; i < D.rows(); ++i)
//     {
//         if (D(i, i) < 0) 
//         {
//             throw std::runtime_error("Matrix contains negative eigenvalues, cannot compute square root");
//         }
//         D(i, i) = std::sqrt(D(i, i));
//     }

//     return V * D * V.transpose();
// }

// TEST_CASE("c2")
// {

//  Eigen::MatrixXd p(9, 9);
//     p << 9999.94,    0,       0,    17.134,    0,       0,   17.1341,    0,       0,
//           0,      40000,     0,       0,       0,       0,       0,       0,       0,
//           0,         0,      0,       0,       0,       0,       0,       0,       0,
//       17.134,       0,       0,    182.906,    0,       0,   0.0299046,  0,       0,
//           0,         0,      0,       0,    40000,     0,       0,       0,       0,
//           0,         0,      0,       0,       0,       0,       0,       0,       0,
//       17.1341,      0,       0, 0.0299046,    0,       0,    182.907,    0,       0,
//           0,         0,      0,       0,       0,       0,       0,    40000,     0,
//           0,         0,      0,       0,       0,       0,       0,       0,       0;

//     try {
//         Eigen::MatrixXd res = sqrtMatSpectral(p);
//         std::cout << "Square root of matrix:\n" << res << std::endl;
//         std::cout << "multiply matrix:\n" << res*res << std::endl;
//     } catch (const std::exception& e) {
//         std::cout << e.what() << std::endl;
//     }

//     BENCHMARK("bench2"){

//     };
// }