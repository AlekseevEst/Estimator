// #include <iostream>
// #include <Eigen/Dense>
// #include <vector>
// #include <catch2/catch.hpp>
// #include "utils.h"

// using namespace Catch::Benchmark;
// using namespace Eigen;

// #define BIG 1.79769e+308 // max value for double

// constexpr double EPS_D = 1.0e-8; ///< Абсолютная точность по умолчанию при сравнениях чисел типа double (1.0e-8)

// struct LAPJVCsp
// {
//     enum TSearchParam
//     {
//         SP_Min, ///< Поиск минимума
//         SP_Max  ///< Поиск максимума
//     };
//     struct CMatrixCSR
//     {
//         std::vector<double> csr_val; ///< Вектор ненулевых элементов матрицы A[n,m] (n - число строк, m - число столбцов), размер nnz
//         std::vector<int> csr_kk;     ///< Вектор индексов колонок ненулевых элементов, размер равен количеству ненулевых элементов nnz
//         std::vector<int> csr_first;  ///< Вектор начальных смещений в векторе CSR, размер n+1
//     };



//     inline bool IsZeroAbs( double value, const double &eps = EPS_D )
// {
//     return ( std::abs( value ) <= eps );
// }

//     void MatrixDenseToCSC( const MatrixXd &A, std::vector<double> &csc_val, std::vector<int> &csc_kk,
//     std::vector<int> &csc_first )
// {
//     csc_val.clear();
//     csc_kk.clear();
//     csc_first.clear();

//     int nnz = A.nonZeros(); // Число ненулевых элементов
//     csc_val.reserve( nnz );
//     csc_kk.reserve( nnz );
//     csc_first.reserve( A.cols() + 1 );

//     int nnz_in_col = 0; // Кол-во ненулевых элементов (non-zero) в столбце
//     csc_first.push_back(0);// Первый элемент надо занулить

//     for( unsigned long long j = 0; j < A.cols(); j++ ) { // Cтолбцы
//         nnz_in_col = 0;
//         for( unsigned long long i = 0; i < A.rows(); i++ ) { // Cтроки
//             if( !IsZeroAbs( A(i,j) ) ) { // if( A[i,j] != 0 )
//                 csc_val.push_back( A(i,j) );//CSC[nnz] = A(i,j);
//                 csc_kk.push_back( i );
//                 nnz_in_col++;
//             }
//         }
//         csc_first.push_back( csc_first[j] + nnz_in_col );
//     }
// }

//     void updateDual(int nc, VectorXd &d, VectorXd &v, VectorXi &todo, int last, double min_)
//     {
//         for (int k = last; k < nc; k++)
//         {
//             int j0 = todo(k);
//             v(j0) += (d(j0) - min_);
//         }
//     }

//     //----------------------------------------------------------------------------------------------------------------------
//     void updateAssignments(VectorXi &lab, VectorXi &y, VectorXi &x, int j, int i0)
//     {
//         int tmp;
//         while (true)
//         {
//             int i = lab(j);
//             y(j) = i;
//             //(j, x[i]) = (x[i], j);
//             tmp = j;
//             j = x[i];
//             x[i] = tmp;
//             if (i == i0)
//             {
//                 return;
//             }
//         }
//     }

//     //----------------------------------------------------------------------------------------------------------------------
//     int solveForOneL(std::vector<double> &cc_, const std::vector<int> &kk, const std::vector<int> &first,
//                      int l, int nc, VectorXd &d, VectorXi &ok, VectorXi &free, VectorXd &v, VectorXi &lab, VectorXi &todo,
//                      VectorXi &y, VectorXi &x, int td1, double resolution, double infValue, bool &fail)
//     {
//         for (int jp = 0; jp < nc; jp++)
//         {
//             d(jp) = infValue;
//             ok(jp) = 0; // false
//         }
//         double min_ = infValue;
//         int i0 = free(l);
//         int j;
//         for (int t = first[i0]; t < first[i0 + 1]; t++)
//         {
//             j = kk[t];
//             double dj = cc_[t] - v(j);
//             d(j) = dj;
//             lab(j) = i0;
//             //        if( dj <= min_ ) { //POSSIBLE FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//             //        if( ( dj < min_ ) || ( std::abs( dj - min_ ) < resolution ) ) { //POSSIBLE FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//             if (((min_ - dj) > resolution) || (std::abs(dj - min_) < resolution))
//             { // POSSIBLE FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//                 // if( dj < min_ ) {
//                 if ((min_ - dj) > resolution)
//                 {
//                     td1 = -1;
//                     min_ = dj;
//                 }
//                 todo(++td1) = j;
//             }
//         }
//         for (int hp = 0; hp <= td1; hp++)
//         {
//             j = todo(hp);
//             if (y(j) == -1)
//             {
//                 updateAssignments(lab, y, x, j, i0);
//                 return td1;
//             }
//             ok(j) = 1; // true
//         }
//         int td2 = (nc - 1);
//         int last = nc;
//         while (true)
//         {
//             if (td1 < 0)
//             {
//                 fail = true; // FAIL!!!
//                 return 1;
//             }
//             int j0 = todo(td1--);
//             int i = y(j0);
//             todo(td2--) = j0;
//             int tp = first[i];
//             while (kk[tp] != j0)
//             {
//                 tp++;
//             }
//             double h = cc_[tp] - v(j0) - min_;
//             for (int t = first[i]; t < first[i + 1]; t++)
//             {
//                 j = kk[t];
//                 //            if( !ok(j) ) {
//                 if (ok(j) == 0)
//                 { // if( false )
//                     double vj = cc_[t] - v(j) - h;
//                     //                if( vj < d(j) ) { // POSSIBLE FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//                     if ((d(j) - vj) > resolution)
//                     { // POSSIBLE FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//                         d(j) = vj;
//                         lab(j) = i;
//                         //                    if( vj == min_ ) { // POSSIBLE FLOWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
//                         if (std::abs(vj - min_) < resolution)
//                         {
//                             if (y[j] == -1)
//                             {
//                                 updateDual(nc, d, v, todo, last, min_);
//                                 updateAssignments(lab, y, x, j, i0);
//                                 return td1;
//                             }
//                             todo(++td1) = j;
//                             ok(j) = 1; // true
//                         }
//                     }
//                 }
//             }
//             if (td1 == -1)
//             {
//                 // The original Pascal code uses finite numbers instead of double.PositiveInfinity
//                 // so we need to adjust slightly here.
//                 min_ = infValue;
//                 last = td2 + 1;
//                 for (int jp = 0; jp < nc; jp++)
//                 {
//                     //                if( ( ( d[jp] < min_ ) || ( std::abs( d[jp] - min_ ) < resolution ) ) && !ok(jp) ) {
//                     //                if( ( ( d[jp] < min_ ) || ( std::abs( d[jp] - min_ ) < resolution ) ) && ( ok(jp) == 0 ) ) {
//                     if ((std::abs(d[jp] - infValue) > resolution) && (((min_ - d[jp]) > resolution) || (std::abs(d[jp] - min_) < resolution)) && (ok(jp) == 0))
//                     {
//                         // if( d[jp] < min_ ) {
//                         if ((min_ - d[jp]) > resolution)
//                         {
//                             td1 = -1;
//                             min_ = d(jp);
//                         }
//                         todo(++td1) = jp;
//                     }
//                 }
//                 for (int hp = 0; hp <= td1; hp++)
//                 {
//                     j = todo(hp);
//                     if (y(j) == -1)
//                     {
//                         updateDual(nc, d, v, todo, last, min_);
//                         updateAssignments(lab, y, x, j, i0);
//                         return td1;
//                     }
//                     ok(j) = 1; // true;
//                 }
//             }
//         }
//     }

//     int JVCsparse(const std::vector<double> &cc, const std::vector<int> &kk, const std::vector<int> &first,
//                       TSearchParam sp, double infValue, double resolution, MatrixXd &rowsol, double &lapcost)
//     {
//         // Объявления
//         int nr = first.size() - 1; // Кол-во строк
//         int max_kk = -1;
//         for (unsigned i = 0; i < kk.size(); i++)
//         { // Поиск максимального элемента в массиве kk - это и будет кол-во столбцов
//             if (kk[i] > max_kk)
//             {
//                 max_kk = kk[i];
//             }
//         }
//         int nc = max_kk + 1; // Кол-во столбцов
//         VectorXi x = Eigen::VectorXi::Zero(nr);
//         VectorXi y = Eigen::VectorXi::Zero (nc);
//         VectorXd u = Eigen::VectorXd::Zero(nr);
//         VectorXd v = Eigen::VectorXd::Zero(nc);
//         VectorXd d = Eigen::VectorXd::Zero(nc);
//         VectorXi ok = Eigen::VectorXi::Zero(nc);  
//         VectorXi xinv = Eigen::VectorXi::Zero(nr);
//         VectorXi free = Eigen::VectorXi::Zero(nr);
//         VectorXi todo = Eigen::VectorXi::Zero(nc);
//         VectorXi lab = Eigen::VectorXi::Zero(nc);
  
//         int l0 = 0;

//         x.array() -= 1;
//         y.array() -= 1;
//         free.array() -= 1;
//         todo.array() -= 1;

//         // Поиск минимума/максимума
//         std::vector<double> cc_ = cc;
//         if (sp == TSearchParam::SP_Max)
//         {
//             std::transform(cc_.begin(), cc_.end(), cc_.begin(),
//                            std::bind(std::multiplies<double>(), std::placeholders::_1, -1.0)); // Умножим на -1 для поиска максимума
//         }

//         // The initialization steps of LAPJVsp only make sense for square matrices
//         if (nr == nc)
//         {
//             for (int jp = 0; jp < nc; jp++)
//             {
//                 v(jp) = infValue;
//             }
//             for (int i = 0; i < nr; i++)
//             {
//                 for (int t = first[i]; t < first[i + 1]; t++)
//                 {
//                     int jp = kk[t];
//                     //                if( cc_[t] < v(jp) ) {
//                     if ((v(jp) - cc_[t]) > resolution)
//                     {
//                         v(jp) = cc_[t];
//                         y(jp) = i;
//                     }
//                 }
//             }
//             for (int jp = (nc - 1); jp >= 0; jp--)
//             {
//                 int i = y(jp);
//                 if (x(i) == -1)
//                 {
//                     x(i) = jp;
//                 }
//                 else
//                 {
//                     y(jp) = -1;
//                     // Here, the original Pascal code simply inverts the sign of x; as that
//                     // doesn't play too well with zero-indexing, we explicitly keep track of
//                     // uniqueness instead.
//                     xinv(i) = 1;
//                 }
//             }
//             int lp = 0;
//             for (int i = 0; i < nr; i++)
//             {
//                 if (xinv(i))
//                 {
//                     continue;
//                 }
//                 if (x(i) != -1)
//                 {
//                     double min_ = infValue;
//                     int j1 = x(i);
//                     for (int t = first[i]; t < first[i + 1]; t++)
//                     {
//                         int jp = kk[t];
//                         if (jp != j1)
//                         {
//                             //                        if( ( cc_[t] - v(jp) ) < min_ ) {
//                             if ((min_ - (cc_[t] - v(jp))) > resolution)
//                             {
//                                 min_ = (cc_[t] - v(jp));
//                             }
//                         }
//                     }
//                     u(i) = min_;
//                     int tp = first[i];
//                     while (kk[tp] != j1)
//                     {
//                         tp++;
//                     }
//                     v(j1) = cc_[tp] - min_;
//                 }
//                 else
//                 {
//                     free(lp++) = i;
//                 }
//             }
//             for (int tel = 0; tel < 2; tel++)
//             {
//                 int h = 0;
//                 int l0p = lp;
//                 lp = 0;
//                 while (h < l0p)
//                 {
//                     // Note: In the original Pascal code, the indices of the lowest
//                     // and second-lowest reduced costs are never reset. This can
//                     // cause issues for infeasible problems; see https://stackoverflow.com/q/62875232/5085211
//                     int i = free(h++);

//                     //------------------------------------------------------------------------------------------------------
//                     // ORIGINAL SEARCH OF MIN AND SUBMIN
//                     //------------------------------------------------------------------------------------------------------

//                     int j0p = -1;         // Index of minimum
//                     int j1p = -1;         // Index of subminimum
//                     double v0 = infValue; // Value of minimum
//                     double vj = infValue; // Value of subminimum
//                     for (int t = first[i]; t < first[i + 1]; t++)
//                     {
//                         int jp = kk[t];
//                         double dj = cc_[t] - v(jp);
//                         //                    if( dj < vj ) {
//                         if ((vj - dj) > resolution)
//                         {
//                             //                        if( dj >= v0 ) { // POSSIBLE FLOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
//                             //                        if( ( dj > v0 ) || ( std::abs( dj - v0 ) < resolution ) ) {
//                             if (((dj - v0) > resolution) || (std::abs(dj - v0) < resolution))
//                             {
//                                 vj = dj;
//                                 j1p = jp;
//                             }
//                             else
//                             {
//                                 vj = v0;
//                                 v0 = dj;
//                                 j1p = j0p;
//                                 j0p = jp;
//                             }
//                         }
//                     }
//                     // If the index of the column with the largest reduced cost has not been
//                     // set, no assignment is possible for this row.
//                     if (j0p < 0)
//                     {
//                         return 1; // No feasible solution!!!
//                     }
//                     int i0 = y(j0p);

//                     //------------------------------------------------------------------------------------------------------
//                     // MY SEARCH OF MIN AND SUBMIN
//                     //------------------------------------------------------------------------------------------------------
//                     //                // find minimum and second minimum reduced cost over columns.
//                     //                int j0p = -1; // Index of minimum
//                     //                int j1p = -1; // Index of subminimum
//                     //                double v0 = infValue; // Value of minimum
//                     //                double vj = infValue; // Value of subminimum
//                     //                arma::vec i_th_row( nc, arma::fill::zeros );
//                     //                i_th_row.fill( infValue );
//                     //                for( int t = first[i]; t < first[i + 1]; t++ ) {
//                     //                    int jp = kk[t];
//                     //                    double dj = cc_[t] - v(jp);
//                     //                    i_th_row(jp) = dj;
//                     //                }
//                     //                j0p = arma::index_min( i_th_row );  // Index of minimum
//                     //                v0 = i_th_row(j0p);                 // Value of minimum
//                     //                i_th_row(j0p) = infValue;

//                     //                j1p = arma::index_min( i_th_row );  // Index of subminimum
//                     //                vj = i_th_row(j1p);                 // Value of subminimum

//                     //                int i0 = y(j0p);
//                     //------------------------------------------------------------------------------------------------------
//                     // ORIGINAL
//                     //------------------------------------------------------------------------------------------------------

//                     u(i) = vj;
//                     //                if( v0 < vj ) { // FIX
//                     if ((vj - v0) > resolution)
//                     { // MY
//                         v(j0p) += (v0 - vj);
//                     }
//                     else if (i0 != -1)
//                     {
//                         j0p = j1p;
//                         i0 = y(j0p);
//                     }
//                     x(i) = j0p;
//                     y(j0p) = i;
//                     if (i0 != -1)
//                     {
//                         //                    if( v0 < vj ) { // FIX
//                         if ((vj - v0) > resolution)
//                         { // MY
//                             free(--h) = i0;
//                         }
//                         else
//                         {
//                             free(lp++) = i0;
//                         }
//                     }

//                     //------------------------------------------------------------------------------------------------------
//                     // MY SEARCH OF MIN AND SUBMIN
//                     //------------------------------------------------------------------------------------------------------
//                     //                u(i) = vj;
//                     //                if( ( vj - v0 ) > resolution ) { // if( v0 < vj )
//                     //                    // change the reduction of the minimum column to increase the minimum
//                     //                    // reduced cost in the row to the subminimum.
//                     //                    v(j0p) += ( v0 - vj );
//                     //                } else {
//                     //                    if( i0 > -1 ) { // minimum and subminimum equal.
//                     //                        // minimum column j1 is assigned.
//                     //                        // swap columns j1 and j2, as j2 may be unassigned.
//                     //                        j0p = j1p;
//                     //                        i0 = y(j0p);
//                     //                    }
//                     //                }
//                     //                // (re-)assign i to j1, possibly de-assigning an i0.
//                     //                x(i) = j0p;
//                     //                y(j0p) = i;
//                     //                if( i0 > -1 ) {
//                     //                    if( ( vj - v0 ) > resolution ) { // FIX
//                     //                        free(--h) = i0;
//                     //                    } else {
//                     //                        free(lp++) = i0;
//                     //                    }
//                     //                }
//                     //------------------------------------------------------------------------------------------------------
//                 }
//             } // end for( int tel = 0; tel < 2; tel++ )
//             l0 = lp;
//         }
//         else
//         { // end if( nr == nc )
//             l0 = nr;
//             for (int i = 0; i < nr; i++)
//             {
//                 free(i) = i;
//             }
//         }
//         int td1 = -1;
//         for (int l = 0; l < l0; l++)
//         {
//             bool fail = false;
//             td1 = solveForOneL(cc_, kk, first, l, nc, d, ok, free, v, lab, todo, y, x, td1, resolution, infValue, fail);
//             if (fail)
//             {
//                 return 1;
//             }
//         }
//         // Prapare output - rowsol and lapcost.
//         lapcost = 0.0;
//         for (int i = 0; i < nr; i++)
//         { // i - row index
//             rowsol(0,i) = x(i);
//             const int j_ = rowsol(0,i); // j - col index
//             const int start = first[i];
//             const int end = first[i + 1];
//             for (int j = start; j < end; j++)
//             {
//                 if (j_ == kk[j])
//                 {
//                     lapcost += cc[j];
//                     break;
//                 }
//             }
//         }
//         return 0;
//     }

//     int operator()(const MatrixXd &costmatrix, TSearchParam sp, double infValue, double resolution, MatrixXd &rowsol,
//                   double &lapcost)
//     {
//         CMatrixCSR csr;
//         MatrixDenseToCSC(costmatrix, csr.csr_val, csr.csr_kk, csr.csr_first);
//         for (size_t i = 0; i < csr.csr_val.size(); i++)
//         {
//             std::cout<< csr.csr_val[i];
//             std::cout<<" ";

//         }
//             std::cout<< std::endl;
       
//         for (size_t i = 0; i < csr.csr_kk.size(); i++)
//         {
//             std::cout<< csr.csr_kk[i];
//             std::cout<<" ";

//         }
//             std::cout<< std::endl;
//         for (size_t i = 0; i < csr.csr_first.size(); i++)
//         {
//             std::cout<< csr.csr_first[i];
//             std::cout<<" ";

//         }
//             std::cout<< std::endl;

         
//         int result = JVCsparse(csr.csr_val, csr.csr_kk, csr.csr_first, sp, infValue, resolution, rowsol, lapcost);
//         return result;
//     }
// };

// TEST_CASE("jvcsp")
// {

//     // Eigen::MatrixXd costMatrix(3, 3);
//     // costMatrix <<   0, 4, 3,
//     //                 1, 5, 4,
//     //                 3, 2, 2;

//     Eigen::MatrixXd costMatrix(5, 5);
//     costMatrix <<   2, 4, 3, 0, 0,
//                     1, 5, 4, 0, 0,
//                     3, 2, 2, 0, 0,
//                     5, 1, 8, 0, 0,
//                     1, 2, 4, 0, 0;  

//     // Eigen::MatrixXd costMatrix(8, 8);
//     // costMatrix <<   93, 1e6, 91, 1e6, 1e6, 1e5, 1e6, 1e6,
//     //                 1e6, 93, 90, 1e6, 98, 1e6, 1e5,1e6,
//     //                 1e6, 1e6, 91, 90, 1e6 , 1e6, 1e6, 1e5,
//     //                 1e6, 1e6, 1e6,   1e6, 1e6,  1e6, 1e6, 1e6,
//     //                 1e6, 1e6, 1e6,   1e6, 1e6,  1e6, 1e6, 1e6,
//     //                 1e6, 1e6, 1e6,   1e6, 1e6,  1e6, 1e6, 1e6,
//     //                 1e6, 1e6, 1e6,   1e6, 1e6,  1e6, 1e6, 1e6,
//     //                 1e6, 1e6, 1e6,   1e6, 1e6,  1e6, 1e6, 1e6;


//     Eigen::MatrixXd rowsol = Eigen::MatrixXd::Zero(1,3);

//     double lapcost;
//     LAPJVCsp lap;
//     lap(costMatrix, LAPJVCsp::TSearchParam::SP_Min,1e7,1e-7,rowsol,lapcost);
//     std::cout<<"lapcost =";
//     std::cout<<lapcost<<std::endl;
//     std::cout<<"rowsol = ";
//     std::cout<<rowsol(0,0);
//     std::cout<<rowsol(0,1);
//     std::cout<<rowsol(0,2);
//     //  std::cout<<rowsol(0,3);
//     //   std::cout<<rowsol(0,4);
    

//     BENCHMARK("jvc_bench"){
//         // cost rez = lap (assigncost,rowsol, colsol, u, v);
//     };
// }