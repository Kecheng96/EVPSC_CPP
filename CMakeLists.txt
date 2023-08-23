cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_CXX_COMPILER "C:/mingw64/bin/gcc.exe")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
#set(CMAKE_CXX_FLAGS_DEBUG "-municode")

project(EVPSC_CPP)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORTCOMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Debug)

set( EIGEN3_INCLUDE_DIR "C:/mingw64/include/eigen3")
#find_package( Eigen3 3.3 REQUIRED)

include_directories( ${EIGEN3_INCLUDE_DIR} )

add_executable(EVPSC_CPP EVPSC_CPP.cpp Input.h Toolbox.cpp Processes.cpp Polycrystals.cpp Grains.cpp Modes.cpp Systems_Voce.cpp)

