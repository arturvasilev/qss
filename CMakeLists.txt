cmake_minimum_required(VERSION 2.8)

project(qss)

set(CMAKE_CXX_COMPILER "/opt/intel/bin/icc")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qopenmp -std=c++11")
# set(SINK_LD_LIBRARY_PATH /opt/intel/lib/mic)

set(SOURCES graph.cpp)
# include_directories(/usr/include)

add_executable(sifter sifter.cpp ${SOURCES})
add_executable(optimizer optimizer.cpp ${SOURCES})
target_link_libraries(optimizer nlopt_cxx m)
