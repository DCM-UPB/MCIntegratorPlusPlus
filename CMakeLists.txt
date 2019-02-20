cmake_minimum_required (VERSION 3.5)
include(FindPackageHandleStandardArgs)

project (mci LANGUAGES CXX VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USER_CXX_FLAGS}")

if (USE_COVERAGE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

# find packages

find_package(MPI)
if (MPI_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI=1")
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  message(STATUS "MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
  message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
endif()

message(STATUS "Configured CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}") 
message(STATUS "Configured CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# set header / library paths
include_directories(include/) # headers

enable_testing()

# continue with subdirectories
add_subdirectory(src) 
add_subdirectory(test)
add_subdirectory(examples)
