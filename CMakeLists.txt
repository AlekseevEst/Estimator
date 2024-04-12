cmake_minimum_required(VERSION 3.14)
project(estimator)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 


find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(3RD_PARTY_RELEASES_DIR ${CMAKE_SOURCE_DIR}/3rd_party/releases)

include(CTest)
include(3rd_party/Catch2)
include(3rd_party/eigen)



option(CLANG_TIDY_RUN "Run clang-tidy alongside project build" OFF)
find_program(CLANG_TIDY_PROGRAM clang-tidy)

if(CLANG_TIDY_PROGRAM AND CLANG_TIDY_RUN)
    message(found CLANG)
    set(CLANG_TIDY_CHECKS "cppcoreguidelines-*,-cppcoreguidelines-avoid-magic-numbers,clang-analyzer-*,bugprone-*")
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_PROGRAM};-checks=${CLANG_TIDY_CHECKS}")
else()
    message(CLANG not found)
endif()
set(ARCH "native" CACHE STRING "Architecture to optimize for (-march)")
add_subdirectory(src)
add_subdirectory(test)

find_package(pybind11 REQUIRED)
if (pybind11_FOUND)
    message(STATUS "-----Pybind11 found------")
    add_subdirectory(bindpy)
else(pybind11_FOUND)
    message(STATUS "-----Pybind11 not found------")
endif(pybind11_FOUND)
target_compile_options(estimator PRIVATE -Wall -Wextra -pedantic)