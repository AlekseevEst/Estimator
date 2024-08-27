// #include <catch2/catch.hpp>
// #include <iostream>
// #include <Eigen/Dense>
// using namespace Catch::Benchmark;
// using M = Eigen::MatrixXd;

// M computeMixingProbability(const M &p_ij, const M &mu_i)
// {
//     M cj(1, mu_i.cols());

//     for (size_t j = 0; j < mu_i.cols(); j++)
//     {
//         double c = 0.;
//         for (size_t i = 0; i < p_ij.cols(); i++)
//         {
//             c = c + (p_ij(i, j) * mu_i(0, i));
//         }
//         cj(0, j) = c;
//     }

//     M m_ij(p_ij.rows(), p_ij.cols());

//     for (size_t j = 0; j < p_ij.cols(); j++)
//     {
//         for (size_t i = 0; i < p_ij.rows(); i++)
//         {
//             m_ij(i, j) = p_ij(i, j) * mu_i(0, i) / cj(0, j);
//         }
//     }
//     return m_ij;
// }

// TEST_CASE("mixx")
// {

//     M p(3, 3);
//     p << 0.7, 0.2,0.1,
//          0.4, 0.4, 0.2,
//          0.3, 0.3, 0.4;
//     M m(1,3);
//     m << 0.2,0.5,0.3;
//     std::cout << p << std::endl;
//     std::cout << m << std::endl;
//     M mx = computeMixingProbability(p, m);
//     std::cout << mx << std::endl;

//     //Expected Mix
//     // 0.325581  0.121212 0.0833333
//     // 0.465116  0.606061  0.416667
//     // 0.209302  0.272727       0.5

//     BENCHMARK("mixx_bench"){

//     };
// }