cmake_minimum_required(VERSION 3.17)
project(mathlib)

set(CMAKE_CXX_STANDARD 20)

add_executable(mathlib main.cpp)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -march=skylake -ffast-math -static-libstdc++ -DX86_64")
set(CMAKE_CXX_FLAGS_RELEASE "-O0 -march=skylake -ffast-math -static-libstdc++ -DX86_64")
