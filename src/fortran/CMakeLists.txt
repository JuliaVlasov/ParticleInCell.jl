set(default_build_type "Release")
cmake_minimum_required(VERSION 3.9)
enable_language(Fortran)
project( pic DESCRIPTION "fortran kernels for pic")
SET(CMAKE_Fortran_FLAGS_DEBUG "-std=f2008 -fstack-arrays -O0 -g -fbacktrace -Werror=intrinsics-std -Wall \
  -pedantic -Wuninitialized \
  -fcheck=array-temps,bounds,do,pointer,recursion \
  -ffpe-trap=invalid,zero,overflow -Werror")
add_library(pic SHARED pic.f90 )
install(TARGETS pic DESTINATION ${CMAKE_SOURCE_DIR})
