project( pic DESCRIPTION "fortran kernels for pic")
cmake_minimum_required(VERSION 3.9)
set(default_build_type "Release")
enable_language(Fortran)
SET(CMAKE_Fortran_FLAGS_DEBUG "-fstack-arrays -O0 -g -fbacktrace -Werror=intrinsics-std -Wall \
  -pedantic -Wuninitialized \
  -fcheck=array-temps,bounds,do,pointer,recursion \
  -ffpe-trap=invalid,zero,overflow -Werror")
add_library(pic SHARED pic.f90 maxwell.f90 particules.f90 quietstart.f90 zone.f90)
install(TARGETS pic DESTINATION ${CMAKE_SOURCE_DIR})
add_executable( vm2d.exe main.f90)
target_link_libraries( vm2d.exe pic)
