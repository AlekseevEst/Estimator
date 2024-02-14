set(LIBRARY_OUTPUT_DIRECTORY /home/dmitry/work/Estimator/scripts)
pybind11_add_module(${PROJECT_NAME}
    MODULE
        bind_ukf.cpp
        bindings.cpp
        )


target_link_libraries(${PROJECT_NAME} PRIVATE eigen filters)
set_target_properties(${PROJECT_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_DIRECTORY})
