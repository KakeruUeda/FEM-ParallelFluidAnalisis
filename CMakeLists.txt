cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_COMPILER "/opt/petsc/bin/mpic++")
set(base_compile_options ${base_compile_options} " -fopenmp")
#set(CMAKE_CXX_FLAGS_RELEASE "-pg")
string(REPLACE " " " " compile_options ${base_compile_options})

project(cfdMPI C CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)



##################################
# include directories
##################################
include_directories(
  /opt/petsc/include
  /opt/homebrew/include/eigen3
)

##################################
# lib directories
##################################
link_directories(
  /opt/petsc/lib
)


add_subdirectory(src)
