cmake_minimum_required(VERSION 3.10)
project(BAB)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CPP_COMPILER "clang++")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(include)
add_executable(test1
    example/main.cpp
    src/Solve.cpp
    src/simplex.cpp
)
