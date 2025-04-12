set(BASE_TARGET "Base")
file(GLOB_RECURSE BASE_HEADERS *.h)
file(GLOB_RECURSE BASE_SOURCES *.cpp)

add_library(${BASE_TARGET} SHARED ${BASE_SOURCES})
target_include_directories(${BASE_TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${BASE_TARGET} geogram Eigen3::Eigen)