// #include <iostream>
// #include <Eigen/Dense>
// #include <vector>
// #include <limits>
// #include <catch2/catch.hpp>
// using namespace Catch::Benchmark;

// using namespace Eigen;

// #include <iostream>
// #include <Eigen/Dense>
// #include <vector>
// #include <limits>

// using namespace Eigen;
// using namespace std;

// typedef MatrixXd CostMatrix;
// typedef VectorXi RowSolution;
// typedef VectorXi ColSolution;
// typedef VectorXd DualVars;

// const double BIG = std::numeric_limits<double>::max();  // Бесконечность

// // Функция для решения задачи о назначениях
// double lapjv(const CostMatrix &cost, RowSolution &rowsol, ColSolution &colsol, DualVars &u, DualVars &v) {
//     int dim = cost.rows();  // Размер задачи
//     vector<int> freeRows;   // Свободные строки
//     vector<int> pred(dim);  // Предшествующие строки
//     VectorXd d(dim);        // Стоимость пути

//     rowsol = VectorXi::Constant(dim, -1);
//     colsol = VectorXi::Constant(dim, -1);
//     u = VectorXd::Zero(dim);  // Двойственные переменные для строк
//     v = VectorXd::Zero(dim);  // Двойственные переменные для столбцов

//     // Шаг 1: Сокращение по строкам
//     for (int i = 0; i < dim; ++i) {
//         double minCost = cost.row(i).minCoeff();
//         u[i] = minCost;
//         for (int j = 0; j < dim; ++j) {
//             if (cost(i, j) == minCost && colsol[j] == -1) {
//                 colsol[j] = i;
//                 rowsol[i] = j;
//                 break;
//             }
//         }
//         if (rowsol[i] == -1) {
//             freeRows.push_back(i);  // Строка осталась неназначенной
//         }
//     }

//     // Шаг 2: Обработка свободных строк
//     while (!freeRows.empty()) {
//         int freerow = freeRows.back();
//         freeRows.pop_back();

//         // Инициализация поиска пути
//         d = VectorXd::Constant(dim, BIG);
//         VectorXi collist(dim);
//         collist.setLinSpaced(dim, 0, dim - 1);
//         int last = 0, low = 0, up = 0;
//         bool found = false;
//         int endOfPath = -1;

//         // Поиск увеличивающего пути
//         while (!found) {
//             if (low == up) {  // Если нечего сканировать, ищем минимальное значение в d
//                 double min = BIG;
//                 for (int k = low; k < dim; ++k) {
//                     int j = collist[k];
//                     if (d[j] < min) {
//                         min = d[j];
//                         endOfPath = j;
//                     }
//                 }

//                 // Обновляем цены столбцов
//                 for (int k = 0; k < last + 1; ++k) {
//                     int j = collist[k];
//                     v[j] += d[j] - min;
//                 }

//                 // Переходим к новой свободной строке
//                 for (int k = last + 1; k--;) {
//                     int j = collist[k];
//                     int i = pred[j];
//                     colsol[j] = i;
//                     int temp = rowsol[i];
//                     rowsol[i] = j;
//                     endOfPath = temp;
//                 }
//             }

//             // Сканирование столбцов
//             for (int k = low; k < dim; ++k) {
//                 int j = collist[k];
//                 if (rowsol[j] == -1) {
//                     found = true;
//                     endOfPath = j;
//                     break;
//                 } else {
//                     // Обновляем минимальные стоимости
//                     int i = rowsol[j];
//                     double newDist = cost(i, j) - u[i] - v[j];
//                     if (newDist < d[j]) {
//                         d[j] = newDist;
//                         pred[j] = i;
//                         collist[up++] = j;
//                     }
//                 }
//             }
//         }
//     }

//     // Рассчитываем общую стоимость
//     double lapcost = 0.0;
//     for (int i = 0; i < dim; ++i) {
//         lapcost += cost(i, rowsol[i]);
//     }

//     return lapcost;
// }


// TEST_CASE("jvc") {
//     // Пример использования
//     CostMatrix cost(3, 3);
//     cost << 4, 1, 3,
//             2, 0, 5,
//             3, 2, 2;

//     RowSolution rowsol;
//     ColSolution colsol;
//     DualVars u, v;

//     double totalCost = lapjv(cost, rowsol, colsol, u, v);

//     cout << "Минимальная стоимость: " << totalCost << endl;
//     cout << "Назначения строк столбцам: " << rowsol.transpose() << endl;
//     cout << "Назначения столбцов строкам: " << colsol.transpose() << endl;
//     BENCHMARK("jvc_bench"){

//     };
// }