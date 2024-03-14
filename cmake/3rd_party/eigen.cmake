cmake_minimum_required(VERSION 3.14)
include(FetchContent)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

FetchContent_Declare(
    eigen
    URL ${3RD_PARTY_RELEASES_DIR}/eigen-3.4.0.tar.gz
)

FetchContent_GetProperties(eigen)

if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
    add_library(eigen INTERFACE)
    target_include_directories(eigen SYSTEM INTERFACE ${eigen_SOURCE_DIR})
endif()
