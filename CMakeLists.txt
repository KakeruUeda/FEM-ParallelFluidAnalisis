cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_COMPILER mpic++)
set(base_compile_options ${base_compile_options} "-O2")
#set(CMAKE_CXX_FLAGS_RELEASE "-pg")
string(REPLACE " " " " compile_options ${base_compile_options})

project(cfdMPI C CXX)
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)


##################################
         # Display #
##################################
message(" ")
message("${PROJECT_NAME}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "COMPILE_OPTIONS       : "      ${compile_options})
message(" ")

##################################
     # include directories #
##################################
include_directories(
  /opt/homebrew/Cellar/metis/5.1.0/include
  /opt/homebrew/include/eigen3
  /usr/local/TextParser/include
  /usr/local/petsc/include
)

##################################
      # lib directories #
##################################
link_directories(
  /opt/homebrew/Cellar/metis/5.1.0/lib
  /usr/local/TextParser/lib
  /usr/local/petsc/lib
  /opt/homebrew/opt/libomp/lib
)


add_subdirectory(src)
