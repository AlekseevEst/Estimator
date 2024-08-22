
#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkfImm_ConteinerCVCACTxz
{
private:
    Track<Eigen::MatrixXd, Imm<Eigen::MatrixXd, ConteinerCVCACTxz>,
          InitUkfImmMeasureModelSph<Eigen::MatrixXd, ConteinerCVCACTxz>>
        track;

public:

    BindTrackUkfImm_ConteinerCVCACTxz(const Detection<Eigen::MatrixXd>& detection):track(detection){}

    Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
    {
        return track.step(detection);
    }

    Eigen::MatrixXd step(double dt)
    {
        return track.step(dt);
    }
    
    Eigen::MatrixXd get_m_i()
    {
        return track.estimator->mu_i;
    }
};

void bind_track_Imm_ConteinerCVCACTxz(pybind11::module &m)
{
    py::class_<BindTrackUkfImm_ConteinerCVCACTxz>(m, "BindTrackUkfImm_ConteinerCVCACTxz")
        .def(py::init<const Detection<Eigen::MatrixXd>&>())
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm_ConteinerCVCACTxz::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkfImm_ConteinerCVCACTxz::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm_ConteinerCVCACTxz::*)(double)) & BindTrackUkfImm_ConteinerCVCACTxz::step)
        .def("get_m_i", (Eigen::MatrixXd(BindTrackUkfImm_ConteinerCVCACTxz::*)()) & BindTrackUkfImm_ConteinerCVCACTxz::get_m_i);
}
