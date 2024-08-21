
#include "bind_track.h"
namespace py = pybind11;

class BindTrackUkfImm_CTxz
{
private:
    TrackImm<Eigen::MatrixXd, Imm<Eigen::MatrixXd, ConteinerCVCACTxz>,
          InitUkfImmCVCTxzCAMeasureModelSph<Eigen::MatrixXd, ConteinerCVCACTxz>>
        track;

public:

    BindTrackUkfImm_CTxz(const Detection<Eigen::MatrixXd>& detection):track(detection){}

    Eigen::MatrixXd step(const Detection<Eigen::MatrixXd>& detection)
    {
        return track.step(detection);
    }

    Eigen::MatrixXd step(double dt)
    {
        // return track.step(dt);
    }
    
    Eigen::MatrixXd get_m_i()
    {
        return track.estimator->mu_i;
    }
};

void bind_track_Imm_CTxz(pybind11::module &m)
{
    py::class_<BindTrackUkfImm_CTxz>(m, "BindTrackUkfImm_CTxz")
        .def(py::init<const Detection<Eigen::MatrixXd>&>())
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm_CTxz::*)(const Detection<Eigen::MatrixXd>&)) & BindTrackUkfImm_CTxz::step)
        .def("step", (Eigen::MatrixXd(BindTrackUkfImm_CTxz::*)(double)) & BindTrackUkfImm_CTxz::step)
        .def("get_m_i", (Eigen::MatrixXd(BindTrackUkfImm_CTxz::*)()) & BindTrackUkfImm_CTxz::get_m_i);
}
