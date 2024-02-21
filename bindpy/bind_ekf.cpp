// #include "bind_ekf.h"
// namespace py = pybind11;

// class BindEkf
// {
// private:
//     ExtendedKalmanFilter<> ekf;
// public:
//     Eigen::MatrixXd predEkf(Eigen::MatrixXd &X)
//     {
//         return ekf.predictEkf(X);
//     }
//     Eigen::MatrixXd corrEkf(Eigen::MatrixXd &Z)
//     {
//        return ekf.correctEkf(Z);
//     }
// };

// void bind_ekf(pybind11::module &m)
// {
//     py::class_<BindEkf>(m, "BindEkf")
//         .def(py::init<>())
//         .def("predictEkf" ,&BindEkf::predEkf)
//             // py::arg("X"))
//         .def("correctEkf" ,&BindEkf::corrEkf);
//             // py::arg("Z"));
// }