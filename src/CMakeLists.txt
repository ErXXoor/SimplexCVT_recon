project(src)
set(CMAKE_CXX_STANDARD 14)
file(GLOB_RECURSE SOURCE_FILES *.c *.cpp)
file(GLOB_RECURSE HEADER_FILES *.h *.hpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

#link to base
target_include_directories(${PROJECT_NAME} PRIVATE Base)
target_link_libraries(${PROJECT_NAME} Base)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} geogram Eigen3::Eigen)