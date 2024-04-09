
####Требования к сборке
option(ESTIMATOR_ADD_TESTS "Add estimator tests to CTest" ON)

#### Добавление исполняемого файла ...#####################################
set(sources_list_test test_estimators.cpp test_angleVel_zero.cpp)
add_executable(estimator_unit_tests ${sources_list_test})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
if(Threads_FOUND)
    message("---Threads FOUND---")
else(message("---Threads not FOUND---"))
endif(Threads_FOUND)

#### Связывание библиотек eigen, Catch2 с ${PROJECT_NAME} и тестом ...#####################################
target_link_libraries(estimator_unit_tests PRIVATE Catch2 eigen filters)

#### Добавление тестов ...#####################################
if(ESTIMATOR_ADD_TESTS)
    include(CTest)
    include(Catch)
    catch_discover_tests(estimator_unit_tests)

    find_program(VALGRIND_PROGRAM NAMES valgrind)
    if(VALGRIND_PROGRAM)
        set(VALGRIND_OPTIONS
            --leak-check=full
            --errors-for-leak-kinds=all
            --show-leak-kinds=all
            --error-exitcode=1
        )
        add_test(
            NAME Valgrind_1
            COMMAND ${VALGRIND_PROGRAM} ${VALGRIND_OPTIONS} ./estimator_unit_tests
            )
    endif()
endif()

