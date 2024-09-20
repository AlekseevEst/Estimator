
// #include "bind_track.h"
// namespace py = pybind11;

// class BindTrackUkfImm_CVCTxyCA
// {
// private:
//     Track<Eigen::MatrixXd, IMM<Eigen::MatrixXd, InitImmFilter1<Eigen::MatrixXd>>>
//         track;

// public:
//     void init(const Detection<Eigen::MatrixXd>& detection)
//     {
//         track.Initialization(detection);
//     }

//     Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
//     {
//         return track.step(detection);
//     }

//     Eigen::MatrixXd step(double dt)
//     {
//         return track.step(dt);
//     }
// };

// void bind_track_Imm_ConteinerCVCTxyCA(pybind11::module &m)
// {
//     py::class_<BindTrackUkfImm_CVCTxyCA>(m, "BindTrackUkfImm_CVCTxyCA")
//         .def(py::init<>())
//         .def("init", (void(BindTrackUkfImm_CVCTxyCA::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkfImm_CVCTxyCA::init)
//         .def("step", (Eigen::MatrixXd(BindTrackUkfImm_CVCTxyCA::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkfImm_CVCTxyCA::step)
//         .def("step", (Eigen::MatrixXd(BindTrackUkfImm_CVCTxyCA::*)(double)) & BindTrackUkfImm_CVCTxyCA::step);

// }
