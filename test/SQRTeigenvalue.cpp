// #include <catch2/catch.hpp>
// #include <iostream>
// #include <Eigen/Dense>
// using namespace Catch::Benchmark;
// using M = Eigen::MatrixXd;

// bool isPositiveDefinite(const M &matrix)
// {
//     using namespace Eigen;
//     Eigen::SelfAdjointEigenSolver<M> eigenSolver(matrix);
//     if (eigenSolver.info() != Eigen::Success)
//     {
//         std::cout << "eigenSolver error";
//         return false;
//     }

//     // Check if all the eigenvalues are positive
//     VectorXd eigenvalues = eigenSolver.eigenvalues();

//     for (int i = 0; i < eigenvalues.size(); ++i)
//     {
//         std::cout << eigenvalues[i] << std::endl;
//         if (eigenvalues[i] <= 0)
//         {
//             std::cout << "eigenvalues отрицательное или 0.";
//         }
//     }
//     return true;
// }

// M sqrtMat(const M &P)
// {
//     Eigen::LLT<M> lltofP(P);
//     if (lltofP.info() != Eigen::Success)
//     {
//         throw std::runtime_error("cholesky decomposition ERROR");
//     }
//     M L = lltofP.matrixL();
//     return L;
// }

// M sqrtMatSpectral(const M &P)
// {
//     Eigen::SelfAdjointEigenSolver<M> eigensolver(P);
//     if (eigensolver.info() != Eigen::Success)
//     {
//         throw std::runtime_error("Eigen decomposition ERROR");
//     }

//     M D = eigensolver.eigenvalues().asDiagonal();
//     std::cout<<D<<std::endl;;
//     M V = eigensolver.eigenvectors();

//     for (int i = 0; i < D.rows(); ++i)
//     {
//         // if (D(i, i) < 0)
//         // {
//         //     if (std::abs(D(i, i)) < 1e-11) // Порог для малых отрицательных значений
//         //     {
//         //         D(i, i) = 0.0;
//         //     }
//         //     else
//         //     {
//         //         throw std::runtime_error("Matrix contains significant negative eigenvalues, cannot compute square root");
//         //     }
//         // }
//         D(i, i) = std::sqrt(D(i, i));
//     }

//     return V * D * V.transpose();
// }

// TEST_CASE("sqrt")
// {
//     M p1(6, 6);
//     M p(6, 6);
//     M p2(7, 7);
//     // p << 1794.26, 0., 3268.7, 0., 1629.12, 0.,
//     //     0, 40000., 0., 0., 0., 0.,
//     //     3268.7, 0., 6735.03, 0., 3273.3, 0.,
//     //     0., 0., 0., 40000., 0., 0.,
//     //     1629.12, 0., 3273.3, 0., 1839.17, 0.,
//     //     0., 0., 0., 0., 0., 40000.;

//     // p1 << 1794.26, 0., 3268.7, 0., 1629.12, 0.,
//     //     0, 40000., 0., 0., 0., 0.,
//     //     3268.7, 0., 6735.03, 0., 3273.3, 0.,
//     //     0., 0., 0., 40000., 0., 0.,
//     //     1629.12, 0., 3273.3, 0., 1839.17, 0.,
//     //     0., 0., 0., 0., 0., 40000.;

//     p2 << 1794.26, 0., 3268.7, 0., 1629.12, 0., 0.,
//         0, 40000., 0., 0., 0., 0., 0.,
//         3268.7, 0., 6735.03, 0., 3273.3, 0., 0.,
//         0., 0., 0., 40000., 0., 0., 0.,
//         1629.12, 0., 3273.3, 0., 1839.17, 0., 0.,
//         0., 0., 0., 0., 0., 40000., 0.,
//         0., 0., 0., 0., 0., 0., 0.;

//         isPositiveDefinite(p2);

//     M q1 = sqrtMatSpectral(p2);
//     std::cout << (q1 * q1.transpose()) << std::endl;

//     std::cout<< p2.determinant()<< std::endl; // детерминант равен нулю, если есть нулевые столбцы и строки

//     // M q = sqrtMat(p2); // холецкий возможен только при отсутсвие нулевых строк и столбцов
//     // std::cout << (q * q.transpose()) << std::endl;

//     BENCHMARK("sqrt_bench"){

//     };
// }
