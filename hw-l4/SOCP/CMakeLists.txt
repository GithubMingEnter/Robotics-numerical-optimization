cmake_minimum_required(VERSION 3.24)
project(SOCP)

set(CMAKE_CXX_STANDARD 14)
find_package(Eigen3 REQUIRED)

include_directories(
        include
        ${EIGEN3_INCLUDE_DIRS}
)
add_executable(SOCP main.cpp)
