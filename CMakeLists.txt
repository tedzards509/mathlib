cmake_minimum_required(VERSION 3.17)
project(mathlib)

find_package(Lua REQUIRED)

set(CMAKE_CXX_STANDARD 20)

add_executable(mathlib main.cpp)
add_executable(mathliblua aml_lua_binding.cpp aml_lua_binding.h testlua.cpp)

target_link_libraries(mathliblua lua)


if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

set(CMAKE_C_COMPILER "clang")
set(CMAKE_CXX_COMPILER "clang++")

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG -O3 -march=native -ffast-math")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")
