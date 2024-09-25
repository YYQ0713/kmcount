cmake_minimum_required(VERSION 2.8)
project(kmc)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (CMAKE_VERSION VERSION_LESS "3.1")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else ()
    set(CMAKE_CXX_STANDARD 11)
endif ()

include_directories(include)
FILE(GLOB_RECURSE KMC "src/*.cpp")
add_subdirectory(libbloom)

find_package(ZLIB REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-function -lz -lm")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ZLIB_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O1 -D_LIBCPP_DEBUG -D_GLIBCXX_DEBUG")


add_executable(kmc ${KMC})
target_link_libraries(kmc ${ZLIB_LIBRARIES})
target_link_libraries(kmc PRIVATE LIBBLOOM)