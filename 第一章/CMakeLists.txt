cmake_minimum_required(VERSION 3.0.2)
project(line_search)

add_compile_options(-std=c++17)
# 需要生成release,运行可执行文件时会出现如下错误
#set(CMAKE_BUILD_TYPE "Release")

find_package(Eigen3 REQUIRED)
if((Eigen3_FOUND) AND (DEFINED EIGEN3_INCLUDE_DIR))
 message(STATUS "Found Eigen3: " ${EIGEN3_INCLUDE_DIR})
 INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
else()
 message(FATAL_ERROR "EIGEN3 not found.")
endif()


add_executable(line_search main.cpp)
add_executable(line_search_raw main_raw.cpp)

