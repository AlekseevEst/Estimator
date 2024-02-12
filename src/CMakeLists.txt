# Список всех .cpp файлов в папке src
file(GLOB SOURCES *.cpp *h)

add_library(filters STATIC ${SOURCES})
target_include_directories(filters PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(filters PRIVATE eigen)