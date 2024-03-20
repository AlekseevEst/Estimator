set(LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../scripts)

pybind11_add_module(${PROJECT_NAME}
    MODULE
        bind_ukf.cpp
        bind_ekf.cpp
        bind_track.cpp
        bindings.cpp
        )


target_link_libraries(${PROJECT_NAME} PRIVATE eigen filters)
set_target_properties(${PROJECT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_DIRECTORY})
