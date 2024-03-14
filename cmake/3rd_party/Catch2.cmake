cmake_minimum_required(VERSION 3.14)
include(FetchContent)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

FetchContent_Declare(
    Catch2
    URL ${3RD_PARTY_RELEASES_DIR}/Catch2-2.13.9.tar.gz
)
FetchContent_MakeAvailable(Catch2)

FetchContent_GetProperties(Catch2 SOURCE_DIR CATCH2_SOURCE_DIR)
list(APPEND CMAKE_MODULE_PATH ${CATCH2_SOURCE_DIR}/contrib)
target_compile_definitions(Catch2 INTERFACE CATCH_CONFIG_ENABLE_BENCHMARKING)
