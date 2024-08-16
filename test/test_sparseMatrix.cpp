// #include <catch2/catch.hpp>
// #include <Eigen/Sparse>
// #include <iostream>
// using namespace Catch::Benchmark;


// TEST_CASE("SparseMatrix")
// {

//         typedef Eigen::SparseMatrix<int> SpMat;
//         typedef Eigen::Triplet<int> T;


//         SpMat Hp(3,5);
//         std::vector<T> tripletList;
//         tripletList.reserve(3);

//         tripletList.push_back(T(0, 0, 1));
//         tripletList.push_back(T(1, 2, 1));
//         tripletList.push_back(T(2, 4, 1));
//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//         Eigen::MatrixXi A (3,1);

//         A << 100,
//              100,
//              100;

//         // std::cout<<Hp<<std::endl<<std::endl;
//         Eigen::MatrixXi X0 = Hp.transpose() * A;

//         Eigen::MatrixXi B (5,1);
//         B << 100,
//               0,
//              100,
//               0,
//              100;
//         // std::cout<<X0<<std::endl<<std::endl;
//         // std::cout<<B<<std::endl<<std::endl;

//         CHECK(X0 == B);
        
//     BENCHMARK("SparseMatrix_multiplication")
//     {
//         Hp.transpose() * A;
//     };
// }

// TEST_CASE("TransforMatState")
// {

//         typedef Eigen::SparseMatrix<int> SpMat;
//         typedef Eigen::Triplet<int> T;


//         SpMat Hp(9,6);
//         std::vector<T> tripletList;
//         tripletList.reserve(6);

//         tripletList.push_back(T(0, 0, 1));
//         tripletList.push_back(T(1, 1, 1));
//         tripletList.push_back(T(3, 2, 1));
//         tripletList.push_back(T(4, 3, 1));
//         tripletList.push_back(T(6, 4, 1));
//         tripletList.push_back(T(7, 5, 1));
              

//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//         Eigen::MatrixXi A (6,1);

//         A << 100,
//              100,
//              100,
//              100,
//              100,
//              100;

//         // std::cout<<Hp<<std::endl<<std::endl;
//         Eigen::MatrixXi X0 = Hp * A;

//         Eigen::MatrixXi B (9,1);
//         B << 100,
//              100,
//               0,
//              100,
//              100,
//               0,
//              100,
//              100,
//               0;
//         // std::cout<<X0<<std::endl<<std::endl;
//         // std::cout<<B<<std::endl<<std::endl;

//         CHECK(X0 == B);
        
//     BENCHMARK("Transform")
//     {

//     };
// }

// TEST_CASE("TransforMatCov")
// {
  

//     // Создаем исходную матрицу 6x6 и заполняем её случайными значениями
//             Eigen::MatrixXi mat6x6 (6,6);

//         mat6x6 << 100, 100, 100, 100, 100, 100,
//              100, 100, 100, 100, 100, 100,
//              100, 100, 100, 100, 100, 100,
//              100, 100, 100, 100, 100, 100,
//              100, 100, 100, 100, 100, 100,
//              100, 100, 100, 100, 100, 100;
//     std::cout << "Original 6x6 Matrix:\n" << mat6x6 << "\n\n";

//     // Создаем матрицу расширения 9x6 для вставки строк и нулевых строк
//     Eigen::MatrixXi expandMat(9, 6);
//     expandMat << 1, 0, 0, 0, 0, 0,
//                  0, 1, 0, 0, 0, 0,
//                  0, 0, 0, 0, 0, 0,
//                  0, 0, 1, 0, 0, 0,
//                  0, 0, 0, 1, 0, 0,
//                  0, 0, 0, 0, 0, 0,
//                  0, 0, 0, 0, 1, 0,
//                  0, 0, 0, 0, 0, 1,
//                  0, 0, 0, 0, 0, 0;

//     // Расширяем исходную матрицу с помощью матрицы расширения
//     Eigen::MatrixXi expandedMat9x6 = expandMat * mat6x6;

//     // Создаем финальную матрицу 9x9, добавляя нулевые столбцы
//     Eigen::MatrixXi mat9x9 = Eigen::MatrixXi::Zero(9, 9);
//     mat9x9.leftCols(6) = expandedMat9x6;

//     std::cout << "Expanded 9x9 Matrix:\n" << mat9x9 << "\n";


//     BENCHMARK("TransformCov")
//     {

//     };
// }


// TEST_CASE("TransforMatStateCVToCA")
// {

//         typedef Eigen::SparseMatrix<double> SpMat;
//         typedef Eigen::Triplet<double> T;


//         SpMat Hp(6,9);
//         std::vector<T> tripletList;
//         tripletList.reserve(6);

//         tripletList.push_back(T(0, 0, 1));
//         tripletList.push_back(T(1, 1, 1));
//         tripletList.push_back(T(2, 3, 1));
//         tripletList.push_back(T(3, 4, 1));
//         tripletList.push_back(T(4, 6, 1));
//         tripletList.push_back(T(5, 7, 1));
              

//         Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//         Eigen::MatrixXd A (6,1);

//         A << 100,
//              100,
//              100,
//              100,
//              100,
//              100;

//         std::cout<<Hp<<std::endl<<std::endl;
//         Eigen::MatrixXd X0 = Hp.transpose() * A;

//         Eigen::MatrixXd B (9,1);
//         B << 100,
//              100,
//               0,
//              100,
//              100,
//               0,
//              100,
//              100,
//               0;
//         std::cout<<X0<<std::endl<<std::endl;
//         // std::cout<<B<<std::endl<<std::endl;

//         // CHECK(X0 == B);
        
//     BENCHMARK("Transform")
//     {

//     };
// }

// TEST_CASE("TransforMatStateCToCA")
// {

//      typedef Eigen::SparseMatrix<double> SpMat;
//      typedef Eigen::Triplet<double> T;

//      SpMat Hp(7, 9);
//      std::vector<T> tripletList;
//      tripletList.reserve(6);

//      tripletList.push_back(T(0, 0, 1));
//      tripletList.push_back(T(1, 1, 1));
//      tripletList.push_back(T(2, 3, 1));
//      tripletList.push_back(T(3, 4, 1));
//      tripletList.push_back(T(4, 6, 1));
//      tripletList.push_back(T(5, 7, 1));

//      Hp.setFromTriplets(tripletList.begin(), tripletList.end());

//      Eigen::MatrixXd A(7, 1);

//      A << 100,
//          100,
//          100,
//          100,
//          100,
//          100,
//          100;

//      std::cout << Hp << std::endl
//                << std::endl;
//      Eigen::MatrixXd X0 = Hp.transpose() * A;

//      Eigen::MatrixXd B(7, 1);
//      B << 100,
//          100,
//          100,
//          100,
//          100,
//          100,
//          0;
//      std::cout << X0 << std::endl
//                << std::endl;
//      // std::cout<<B<<std::endl<<std::endl;

//      // CHECK(X0 == B);

//      BENCHMARK("Transform"){

//      };
// }
