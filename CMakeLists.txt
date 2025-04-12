cmake_minimum_required(VERSION 3.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(SimplexCVT_recon)

set(CMAKE_CXX_STANDARD 11)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#include(stan_math)
#include(libigl)
#igl_include(glfw)

find_package(Eigen3 REQUIRED)
if (EIGEN3_FOUND)
    message("Eigen3 found")
endif ()
include(geogram)

add_subdirectory(Base)
add_subdirectory(src)