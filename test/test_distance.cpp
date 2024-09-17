// #include <catch2/catch.hpp>
// #include <iostream>
// #include <typeinfo>
// #include <typeindex>

// #include <cstddef>

// #include "Eigen/Dense"
// #include "Eigen/Sparse"
// #include <cmath>
// #include <vector>
// #include <iostream>
// #include "utils.h"
// #include "models.h"
// #include "converter.h"

// using namespace Catch::Benchmark;


// double distance(Eigen::MatrixXd &Z, Eigen::MatrixXd &Se, Eigen::MatrixXd &Ze) 
//     {
//         Eigen::MatrixXd v = Z - Ze;
//         PRINTM(v.transpose() * Se.inverse() * v);
//         double distance = (v.transpose() * Se.inverse() * v)(0,0);
//         return distance;

//     }

// TEST_CASE("distance")
// {
//     Eigen::MatrixXd Z(2,1);
//     Z<<10,5;
//     Eigen::MatrixXd Ze(2,1);
//     Ze<< 9,4;
//     Eigen::MatrixXd Se(2,2);
//     Se<<2,0,0,2; 

//     std::cout<<distance(Z,Se,Ze);
// }