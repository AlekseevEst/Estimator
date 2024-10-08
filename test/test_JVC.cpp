// #include <iostream>
// #include <Eigen/Dense>
// #include <vector>
// #include <catch2/catch.hpp>
// #include "utils.h"

// #include <Eigen/Sparse>
// using namespace Catch::Benchmark;
// using namespace Eigen;

// using SpMat = Eigen::SparseMatrix<double>;
// using T = Eigen::Triplet<double>;

// #define BIG 1.79769e+308 // max value for double

// using cost = double;
// using row = int;
// using col = int;

// template <class M> 
// struct LAPJV
// {

//     double operator()(M& matrixcost, M& rowsol, M& colsol, M& u, M& v)
//     {

// // input:
// // assigncost - cost matrix
// // output:
// // rowsol     - column assigned to row in solution
// // colsol     - row assigned to column in solution
// // u          - dual variables, row reduction numbers
// // v          - dual variables, column reduction numbers
    
//     M assigncost = matrixcost;


//     if (matrixcost.cols() > matrixcost.rows())
//     {
//         assigncost.resize(matrixcost.cols(),matrixcost.cols());
//         assigncost.setConstant(BIG-1.0);
//         assigncost.block(0,0, matrixcost.rows(), matrixcost.cols()) = matrixcost;
//     }  
//     if (matrixcost.rows() > matrixcost.cols())
//     {
//         assigncost.resize(matrixcost.rows(),matrixcost.rows());
//         assigncost.setConstant(BIG-1.0);
//         assigncost.block(0,0, matrixcost.rows(), matrixcost.cols()) = matrixcost;
//     }  

//     int dim = assigncost.rows();
//     rowsol.resize(1,dim);
//     rowsol.setZero();

//     bool unassignedfound;
//     row i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *freeunassigned;
//     col j, j1, j2, endofpath, last, low, up, *collist, *matches;
//     cost min, h, umin, usubmin, v2, *d;

//     freeunassigned = new row[dim];     // list of unassigned rows.
//     collist = new col[dim];  // list of columns to be scanned in various ways.
//     matches = new col[dim];  // counts how many times a row could be assigned.
//     d = new cost[dim];       // 'cost-distance' in augmenting path calculation.
//     pred = new row[dim];     // row-predecessor of column in augmenting/alternating path.

//     // init how many times a row will be assigned in the column reduction.
//     for (i = 0; i < dim; i++) matches[i] = 0;

//     // COLUMN REDUCTION
//     for (j = dim; j--;)  // reverse order gives better results.
//     {
//         // find minimum cost over rows.
//         min = assigncost(0,j);
//         imin = 0;
//         for (i = 1; i < dim; i++)
//             if (assigncost(i,j) < min) {
//                 min = assigncost(i,j);
//                 imin = i;
//             }
//         v(0,j) = min;
//         int j1 = rowsol(0,imin); 
//         if (++matches[imin] == 1) {
//             // init assignment if minimum row assigned for first time.
//             rowsol(0,imin) = j;
//             colsol(0,j) = imin;
//         } else if (v(0,j) < v(0,j1)){ //v(0,rowsol(0,imin))
//             int j1 = rowsol(0, imin);
//             rowsol(0, imin) = j;
//             colsol(0,j) = imin;
//             colsol(0,j1) = -1;
//         } else
//             colsol(0,j) = -1;  // row already assigned, column not assigned.
//     }

//     // REDUCTION TRANSFER
//     for (i = 0; i < dim; i++)
//         if (matches[i] == 0)  // fill list of unassigned 'free' rows.
//             freeunassigned[numfree++] = i;
//         else if (matches[i] == 1)  // transfer reduction from rows that are assigned once.
//         {
//             j1 = rowsol(0,i);
//             min = BIG;
//             for (j = 0; j < dim; j++)
//                 if (j != j1)
//                     if (assigncost(i,j) - v(0,j) < min) min = assigncost(i,j) - v(0,j);
//             v(0,j1) -= min;
//         }

//     //   AUGMENTING ROW REDUCTION
//     int loopcnt = 0;  // do-loop to be done twice.
//     do {
//         loopcnt++;

//         //     scan all free rows.
//         //     in some cases, a free row may be replaced with another one to be scanned next.
//         k = 0;
//         prvnumfree = numfree;
//         numfree = 0;  // start list of rows still free after augmenting row reduction.
//         while (k < prvnumfree) {
//             i = freeunassigned[k];
//             k++;

//             //       find minimum and second minimum reduced cost over columns.
//             umin = assigncost(i,0) - v(0,0);
//             j1 = 0;
//             usubmin = BIG;
//             for (j = 1; j < dim; j++) {
//                 h = assigncost(i,j) - v(0,j);
//                 if (h < usubmin)
//                     if (h >= umin) {
//                         usubmin = h;
//                         j2 = j;
//                     } else {
//                         usubmin = umin;
//                         umin = h;
//                         j2 = j1;
//                         j1 = j;
//                     }
//             }

//             i0 = colsol(0,j1);
//             if (umin < usubmin)
//                 //         change the reduction of the minimum column to increase the minimum
//                 //         reduced cost in the row to the subminimum.
//                 v(0,j1) = v(0,j1) - (usubmin - umin);
//             else              // minimum and subminimum equal.
//                 if (i0 > -1)  // minimum column j1 is assigned.
//             {
//                 //           swap columns j1 and j2, as j2 may be unassigned.
//                 j1 = j2;
//                 i0 = colsol(0,j2);
//             }

//             //       (re-)assign i to j1, possibly de-assigning an i0.
//             rowsol(0,i) = j1;
//             colsol(0,j1) = i;

//             if (i0 > -1)  // minimum column j1 assigned earlier.
//                 if (umin < usubmin)
//                     //           put in current k, and go back to that k.
//                     //           continue augmenting path i - j1 with i0.
//                     freeunassigned[--k] = i0;
//                 else
//                     //           no further augmenting reduction possible.
//                     //           store i0 in list of free rows for next phase.
//                     freeunassigned[numfree++] = i0;
//         }
//     } while (loopcnt < 2);  // repeat once.

//     // AUGMENT SOLUTION for each free row.
//     for (f = 0; f < numfree; f++) {
//         freerow = freeunassigned[f];  // start row of augmenting path.

//         // Dijkstra shortest path algorithm.
//         // runs until unassigned column added to shortest path tree.
//         for (j = dim; j--;) {
//             d[j] = assigncost(freerow,j) - v(0,j);
//             pred[j] = freerow;
//             collist[j] = j;  // init column list.
//         }

//         low = 0;  // columns in 0..low-1 are ready, now none.
//         up = 0;   // columns in low..up-1 are to be scanned for current minimum, now none.
//                   // columns in up..dim-1 are to be considered later to find new minimum,
//                   // at this stage the list simply contains all columns
//         unassignedfound = false;
//         do {
//             if (up == low)  // no more columns to be scanned for current minimum.
//             {
//                 last = low - 1;

//                 // scan columns for up..dim-1 to find all indices for which new minimum occurs.
//                 // store these indices between low..up-1 (increasing up).
//                 min = d[collist[up++]];
//                 for (k = up; k < dim; k++) {     
//                     j = collist[k];
//                     h = d[j];
//                     if (h <= min) {
//                         if (h < min)  // new minimum.
//                         {
//                             up = low;  // restart list at index low.
//                             min = h;
//                         }
//                         // new index with same minimum, put on undex up, and extend list.
//                         collist[k] = collist[up];
//                         collist[up++] = j;
//                     }
//                 }
//                 // check if any of the minimum columns happens to be unassigned.
//                 // if so, we have an augmenting path right away.
//                 for (k = low; k < up; k++)
//                     if (colsol(0,collist[k]) < 0) {
//                         endofpath = collist[k];
//                         unassignedfound = true;
//                         break;
//                     }
//             }

//             if (!unassignedfound) {
//                 // update 'distances' between freerow and all unscanned columns, via next scanned
//                 // column.
//                 j1 = collist[low];
//                 low++;
//                 i = colsol(0,j1);
//                 h = assigncost(i,j1) - v(0,j1) - min;

//                 for (k = up; k < dim; k++) {  
//                     j = collist[k];
//                     v2 = assigncost(i,j) - v(0,j) - h;
//                     if (v2 < d[j]) {
//                         pred[j] = i;
//                         if (v2 == min)  // new column found at same minimum value
//                             if (colsol(0,j) < 0) {
//                                 // if unassigned, shortest augmenting path is complete.
//                                 endofpath = j;
//                                 unassignedfound = true;
//                                 break;
//                             }
//                             // else add to list to be scanned right away.
//                             else {
//                                 collist[k] = collist[up];
//                                 collist[up++] = j;
//                             }
//                         d[j] = v2;
//                     }
//                 }
//             }
//         } while (!unassignedfound);

//         // update column prices.
//         for (k = last + 1; k--;) {
//             j1 = collist[k];
//             v(0,j1) += d[j1] - min;
//         }

//         // reset row and column assignments along the alternating path.
//         do {
//             i = pred[endofpath];
//             colsol(0,endofpath) = i;
//             j1 = endofpath;
//             endofpath = rowsol(0,i);
//             rowsol(0,i) = j1;
//         } while (i != freerow);
//     }


//     // calculate optimal cost.
//     cost lapcost = 0;
//     int ii;
//     //  for (i = 0; i < dim; i++)
//     for (i = dim; i--;)
//     {
//         j = rowsol(0, i);
//         ii = colsol(0,i);
//         u(0, i) = assigncost(i, j) - v(0, j);

//         if (assigncost(i, j) == BIG-1.0)
//         {
//             rowsol(0,i) = std::numeric_limits<double>::quiet_NaN();

//         }

//         if (assigncost(ii, i) == BIG-1.0)
//         {
//             colsol(0,i) = std::numeric_limits<double>::quiet_NaN();
//         }
//         if (assigncost(i, j) != BIG-1.0)
//             lapcost += assigncost(i, j);
//     }

//     // free reserved memory.
//     delete[] pred;
//     delete[] freeunassigned;
//     delete[] collist;
//     delete[] matches;
//     delete[] d;
//     return lapcost;

//     }

// };


// bool compareMatricesWithNaN(const Eigen::MatrixXd& mat1, const Eigen::MatrixXd& mat2) {
//     // Сначала проверяем, что размеры матриц совпадают
//     if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols()) {
//         return false;
//     }
//     // Проходим по каждому элементу и сравниваем их
//     for (int i = 0; i < mat1.rows(); ++i) {
//         for (int j = 0; j < mat1.cols(); ++j) {
//             double elem1 = mat1(i, j);
//             double elem2 = mat2(i, j);
//             // Проверяем, если оба значения NaN
//             if (std::isnan(elem1) && std::isnan(elem2)) {
//                 continue; // Если оба NaN, то считаем, что они равны
//             }
//             // Проверяем на равенство (включая числа с плавающей запятой)
//             if (elem1 != elem2) {
//                 return false;
//             }
//         }
//     }
//     return true;
// }

// TEST_CASE("jvc") {


// Eigen::MatrixXd assigncost1(3,5); // строки - трассы, столбцы - отметки.
// assigncost1 <<   0,4,3,1,5,
//                  1,5,4,1,7,
//                  3,2,2,5,9;

// int dim1 = assigncost1.cols();

// Eigen::MatrixXd rowsol1;//(1,dim);
// // rowsol1.setZero();
// Eigen::MatrixXd colsol1(1,dim1);
// colsol1.setZero();
// Eigen::MatrixXd u1(1,dim1);
// u1.setZero();
// Eigen::MatrixXd v1(1,dim1);
// v1.setZero();
// LAPJV<Eigen::MatrixXd> lap1;

// cost rez1 = lap1 (assigncost1,rowsol1, colsol1, u1, v1);
// double exprez1 = 3.0;
// CHECK(rez1 == exprez1);

// Eigen::MatrixXd expectedrowsol1(1,dim1);
// expectedrowsol1 << 0,3,2,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN();


//     CHECK(compareMatricesWithNaN(rowsol1,expectedrowsol1) == true);
//     bool condition = (compareMatricesWithNaN(rowsol1,expectedrowsol1) == true);
//     if(!condition){
//         PRINTM(rowsol1);
//         PRINTM(expectedrowsol1);
//     }

// Eigen::MatrixXd assigncost2(5,3); 
// assigncost2 <<   0,4,3,           
//                 1,5,4,
//                 3,2,2,
//                 1,8,7,
//                 6,9,1;


// int dim2 = assigncost2.rows();

// Eigen::MatrixXd rowsol2;//(1,dim2);
// // rowsol2.setZero();
// Eigen::MatrixXd colsol2(1,dim2);
// colsol2.setZero();
// Eigen::MatrixXd u2(1,dim2);
// u2.setZero();
// Eigen::MatrixXd v2(1,dim2);
// v2.setZero();
// LAPJV<Eigen::MatrixXd> lap2;

// cost rez2 = lap2 (assigncost2,rowsol2, colsol2, u2, v2);


// double exprez2 = 3.;
// Eigen::MatrixXd expectedrowsol2(1,dim2);
// expectedrowsol2 << 0,std::numeric_limits<double>::quiet_NaN(),1,std::numeric_limits<double>::quiet_NaN(),2;

//     CHECK(compareMatricesWithNaN(rowsol2,expectedrowsol2) == true);
//     condition = (compareMatricesWithNaN(rowsol2,expectedrowsol2) == true);
//     if(!condition){
//         PRINTM(rowsol2);
//         PRINTM(expectedrowsol2);
//     }

//     CHECK(rez2 == exprez2);

//     BENCHMARK("jvc_bench"){
// // cost rez = lap (assigncost,rowsol, colsol, u, v);
//     };
// }