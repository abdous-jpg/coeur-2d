cmake_minimum_required(VERSION 3.16)
project(coeur2d C CXX)

set(CMAKE_CXX_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_package(HDF5 REQUIRED COMPONENTS C)
find_package(OpenMP REQUIRED)

add_executable(coeur2d src/main.cpp)
target_include_directories(coeur2d PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(coeur2d PRIVATE ${HDF5_LIBRARIES} OpenMP::OpenMP_CXX)
