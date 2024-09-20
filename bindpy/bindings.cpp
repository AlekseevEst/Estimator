#include "bind_ukf.h"
#include "bind_track.h"

PYBIND11_MODULE(estimator, m) {
    // bind_ukf(m);
    bind_track_CTxy(m);
    // bind_track_CTxz(m);
    bind_track_CV(m);
    bind_track_CA(m);
    // bind_track_Imm_ConteinerCVCTxyCA(m);
    // bind_track_Imm_ConteinerCVCACTxz(m);
}
