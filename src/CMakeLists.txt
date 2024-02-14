# Список всех .cpp файлов в папке src
file(GLOB SOURCES *.cpp)
set (HEADER ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(filters STATIC ${SOURCES} ${HEADER})
target_include_directories(filters PUBLIC ${HEADER})
target_link_libraries(filters PRIVATE eigen)
