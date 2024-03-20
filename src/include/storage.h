#pragma once
#include <vector>
#include "track.h"

struct Storage
{

    void getTrack();
    void putTrack();

    private:
    std::vector<Track> tracks;
    std::vector<Track> usedTracks;
    std::vector<Track> UnusedTracks;

};
