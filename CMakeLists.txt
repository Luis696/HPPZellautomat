cmake_minimum_required(VERSION 3.28)
project(HPPZellautomat)
# Find OpenMP package
find_package(OpenMP REQUIRED)
# Add MPI Package to Project
find_package(MPI REQUIRED)
# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")

set(CMAKE_CXX_STANDARD 14)

add_executable(HPPZellautomat main.cpp)
add_executable(HPPZellautomatV1 mainV1.cpp)
# Add libraries for code completion and compiling
target_link_libraries(HPPZellautomat PUBLIC MPI::MPI_CXX)
target_link_libraries(HPPZellautomatV1 PUBLIC MPI::MPI_CXX)