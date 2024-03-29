#pragma once
#include <vector>
#include <unordered_map>

template <class typeTrack>
struct Storage
{

    // void getTrack();
    // void putTrack();

    private:

    std::vector<typeTrack*> tracks; //сложность  - O(1) 

    std::unordered_map<uint64_t, typeTrack*> tracks; //ассоциативный неупорядоченный контейнер. Сложность O(1). У обычного std::map O(log n)

};
