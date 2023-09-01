cmake_minimum_required(VERSION 3.0.0)

#set(CMAKE_CXX_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")

project(EVPSC_CPP)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Debug)

set( EIGEN3_INCLUDE_DIR "/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
set( OMP_INCLUDE_DIR "/opt/homebrew/Cellar/libomp/14.0.6/include")

include_directories( ${EIGEN3_INCLUDE_DIR} ${OMP_INCLUDE_DIR} )

add_executable(EVPSC EVPSC_CPP.cpp Toolbox.cpp Processes.cpp Polycrystals.cpp Grains.cpp Modes.cpp Systems_Voce.cpp global.cpp)
