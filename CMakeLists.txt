cmake_minimum_required(VERSION 3.17)
project(mathlib)

find_package(Lua REQUIRED)

set(CMAKE_CXX_STANDARD 20)

add_executable(mathlib main.cpp aml_lua_binding.cpp aml_lua_binding.h)

target_link_libraries(mathlib lua)


if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

if (UNIX)

set(CMAKE_C_COMPILER "clang")
set(CMAKE_CXX_COMPILER "clang++")

endif (UNIX)

message(STATUS ${CMAKE_SYSTEM_PROCESSOR})

if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "aarch64")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG -O3 -mcpu=cortex-a72 -ffast-math -static-libstdc++")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mcpu=cortex-a72 -ffast-math -static-libstdc++")


else ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG -O3 -march=skylake -ffast-math -static-libstdc++ -DX86_64")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=skylake -ffast-math -static-libstdc++ -DX86_64")

endif ()