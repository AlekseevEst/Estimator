set (HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(filters INTERFACE)
target_include_directories(filters INTERFACE ${HEADER})
target_link_libraries(filters INTERFACE eigen)
