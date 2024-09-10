
// #include "bind_track.h"
// namespace py = pybind11;

// class BindTrackUkfImm_ConteinerCVCACTxy
// {
// private:
//     Track<Eigen::MatrixXd, Imm<Eigen::MatrixXd, ConteinerCVCACTxy>,
//           InitUkfImmMeasureModelSph<Eigen::MatrixXd, ConteinerCVCACTxy>>
//         track;

// public:

//     BindTrackUkfImm_ConteinerCVCACTxy(const Detection<Eigen::MatrixXd>& detection):track(detection){}

//     Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
//     {
//         return track.step(detection);
//     }

//     Eigen::MatrixXd step(double dt)
//     {
//         return track.step(dt);
//     }
    
//     Eigen::MatrixXd get_m_i()
//     {
//         return track.estimator->mu_i;
//     }
// };

// void bind_track_Imm_ConteinerCVCACTxy(pybind11::module &m)
// {
//     py::class_<BindTrackUkfImm_ConteinerCVCACTxy>(m, "BindTrackUkfImm_ConteinerCVCACTxy")
//         .def(py::init<const Detection<Eigen::MatrixXd>&>())
//         .def("step", (Eigen::MatrixXd(BindTrackUkfImm_ConteinerCVCACTxy::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkfImm_ConteinerCVCACTxy::step)
//         .def("step", (Eigen::MatrixXd(BindTrackUkfImm_ConteinerCVCACTxy::*)(double)) & BindTrackUkfImm_ConteinerCVCACTxy::step)
//         .def("get_m_i", (Eigen::MatrixXd(BindTrackUkfImm_ConteinerCVCACTxy::*)()) & BindTrackUkfImm_ConteinerCVCACTxy::get_m_i);
// }
