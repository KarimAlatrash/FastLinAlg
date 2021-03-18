cmake_minimum_required(VERSION 3.17)
project(FastLinAlg)

set(CMAKE_CXX_STANDARD 11)

add_executable(FastLinAlg
        main.cpp
        Linear_solve.hpp
        linalg.cpp)