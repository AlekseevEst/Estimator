#include <catch2/catch.hpp>
#include <Eigen/Sparse>
#include <iostream>
using namespace Catch::Benchmark;


TEST_CASE("SparseMatrix")
{

        typedef Eigen::SparseMatrix<int> SpMat;
        typedef Eigen::Triplet<int> T;


        SpMat Hp(3,5);
        std::vector<T> tripletList;
        tripletList.reserve(3);

        tripletList.push_back(T(0, 0, 1.0));
        tripletList.push_back(T(1, 2, 1.0));
        tripletList.push_back(T(2, 4, 1.0));
        Hp.setFromTriplets(tripletList.begin(), tripletList.end());

        Eigen::MatrixXi A (5,1);

        A << 100.0,
             100.0,
             100.0, 
             100.0, 
             100.0;

        Eigen::MatrixXi X0 = Hp.transpose() * A;

        Eigen::MatrixXi B (5,1);
        B << 100,
              0,
             100,
              0,
             100;
        CHECK(X0 == B);
        
    BENCHMARK("SparseMatrix_multiplication")
    {
        Hp * A;
    };
}

