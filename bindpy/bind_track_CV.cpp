// #include "bind_track.h"
// namespace py = pybind11;

// // Points Points_CV;
// class BindTrackUkf_CV
// {
// private:
//     Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSphCVCT>> track;

// public:
//     BindTrackUkf_CV(Eigen::MatrixXd state,
//                     Eigen::MatrixXd processNoise,
//                     Eigen::MatrixXd measureNoise,
//                     ParamSigmaPoints paramSigmaPoints) : track(state, processNoise, measureNoise, paramSigmaPoints) {}

//     Eigen::MatrixXd step(double dt, const Eigen::MatrixXd &meas)
//     {
//         return track.Step(dt, meas);
//     }

//     Eigen::MatrixXd step(double dt)
//     {
//         return track.Step(dt);
//     }
// };

// void bind_track_CV(pybind11::module &m)
// {
//     py::class_<BindTrackUkf_CV>(m, "BindTrackUkf_CV")
//         .def(py::init<const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, ParamSigmaPoints>())
//         .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(double, const Eigen::MatrixXd &)) & BindTrackUkf_CV::step)
//         .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(double)) & BindTrackUkf_CV::step);

// }


#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CV
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstVel, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>,
          InitUKFStateModelCVMeasureModelSph<Eigen::MatrixXd, FuncConstVel,
                                             FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZ>>
        track;

public:

    BindTrackUkf_CV(const Detection<Eigen::MatrixXd>& detection):track(detection){}

    Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
    {
        return track.step(detection);
    }

    Eigen::MatrixXd step(double dt)
    {
        return track.step(dt);
    }
};

void bind_track_CV(pybind11::module &m)
{
    py::class_<BindTrackUkf_CV>(m, "BindTrackUkf_CV")
        .def(py::init<const Detection<Eigen::MatrixXd>&>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CV::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CV::*)(double)) & BindTrackUkf_CV::step);
}
