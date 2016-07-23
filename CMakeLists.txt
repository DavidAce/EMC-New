cmake_minimum_required(VERSION 2.8.9)
project(EMC-CLion)
 
#Bring the headers, such as Student.h into the project
include_directories(source)
 
#Can manually add the sources using the set command as follows:
#set(SOURCES src/mainapp.cpp src/Student.cpp)
 
#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "source/*.cpp")
 
add_executable(EMC ${SOURCES})
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
target_compile_options(EMC PUBLIC -std=c++11 -fopenmp -Wall -Wfloat-conversion -O3 -march=native -mno-avx -Wno-deprecated-declarations -DEIGEN_NO_DEBUG -DNDEBUG)

