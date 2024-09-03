
// #include "bind_track.h"
// namespace py = pybind11;

// class BindTrackUkf_CTxz
// {
// private:
//     Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurnXZ, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>,
//           InitUKFStateModelCTMeasureModelSph<Eigen::MatrixXd, FuncConstTurnXZ,
//                                              FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>>
//         track;

// public:

//     BindTrackUkf_CTxz(const Detection<Eigen::MatrixXd>& detection):track(detection){}

//     Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
//     {
//         return track.step(detection);
//     }

//     Eigen::MatrixXd step(double dt)
//     {
//         return track.step(dt);
//     }
// };

// void bind_track_CTxz(pybind11::module &m)
// {
//     py::class_<BindTrackUkf_CTxz>(m, "BindTrackUkf_CTxz")
//         .def(py::init<const Detection<Eigen::MatrixXd>&>())
//         .def("step", (Eigen::MatrixXd(BindTrackUkf_CTxz::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CTxz::step)
//         .def("step", (Eigen::MatrixXd(BindTrackUkf_CTxz::*)(double)) & BindTrackUkf_CTxz::step);
// }
