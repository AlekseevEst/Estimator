// #include <catch2/catch.hpp>
// #include <iostream>
// #include <Eigen/Dense>
// using namespace Catch::Benchmark;
// using M = Eigen::MatrixXd;


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

//     M lambda = eigensolver.eigenvalues().asDiagonal();

//     M V = eigensolver.eigenvectors();
//     for (int i = 0; i < lambda.rows(); ++i)
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

//             if (lambda(i, i) > 0)
//             {
//                 lambda(i, i) = std::sqrt(lambda(i, i));
//             }
//             else
//             {
//                 lambda(i, i) = 0; // Если собственное значение нулевое, оставляем его нулевым
//             }
        
//     }

//     return V * lambda;
// }

// TEST_CASE("U")
// {
  
//     M p(9, 9);
   
//     p << 1779.99,       0,       0, 3257.46 ,      0,       0, 1621.22,       0,       0,
//       0,   40000,       0,       0,       0,       0,       0,       0,       0,
//       0,       0,       0,       0,       0,       0,       0,       0,       0,
// 3257.46,       0,       0, 6750.78,       0,       0, 3276.18,       0,       0,
//       0,       0,       0,       0,   40000,       0,       0,       0,       0,
//       0,       0,       0,       0,       0,       0,       0,       0,       0,
// 1621.22,       0,       0, 3276.18,       0,       0, 1839.02,       0,       0,
//       0,       0,       0,       0,       0,       0,       0,   40000,       0,
//       0,       0,       0,       0,       0,       0,       0,       0,       0;

//       M p1(3, 3);
//     p1 << 4, 1, 2,
//          1, 3, 0,
//          2, 0, 5;


//     M q1 = sqrtMatSpectral(p);
//     std::cout << (q1) << std::endl<< std::endl;
//     std::cout << (q1 * q1.transpose()) << std::endl;

//     M q2 = sqrtMatSpectral(p1);
//     std::cout << (q2) << std::endl<< std::endl;
//     std::cout << (q2 * q2.transpose()) << std::endl;


//     // std::cout<< p.determinant()<< std::endl; // детерминант равен нулю, если есть нулевые столбцы и строки

//     // M q = sqrtMat(p2); // холецкий возможен только при отсутсвие нулевых строк и столбцов
//     // std::cout << (q * q.transpose()) << std::endl;

//     BENCHMARK("U_bench"){

//     };
// }


