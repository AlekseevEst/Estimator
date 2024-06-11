
#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkf_CT
{
private:
    Track<Eigen::MatrixXd, UnscentedKalmanfilter<Eigen::MatrixXd, FuncConstTurn, FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>,
          InitUKFStateModelCTMeasureModelSph<Eigen::MatrixXd, FuncConstTurn,
                                             FuncMeasSphCVCT, FuncControlMatrix_XvXYvYZvZW>>
        track;

public:

    BindTrackUkf_CT(const Detection<Eigen::MatrixXd>& detection):track(detection){}

    Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
    {
        return track.step(detection);
    }

    Eigen::MatrixXd step(double dt)
    {
        return track.step(dt);
    }
};

void bind_track_CT(pybind11::module &m)
{
    py::class_<BindTrackUkf_CT>(m, "BindTrackUkf_CT")
        .def(py::init<const Detection<Eigen::MatrixXd>&>())
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CT::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkf_CT::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkf_CT::*)(double)) & BindTrackUkf_CT::step);
}
