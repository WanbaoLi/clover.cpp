cmake_minimum_required(VERSION 3.18)
project(Bridge213 VERSION 2.1.3 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Generate the upstream configuration headers inside the build tree.
set(generated "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(bridge_version "${PROJECT_VERSION}")
file(CONFIGURE OUTPUT "${generated}/version.h"
  CONTENT "#define BRIDGE_VERSION \"${bridge_version}\"\n")

set(qxs_inline "${generated}/lib_alt_QXS/inline")
file(GLOB inline_headers CONFIGURE_DEPENDS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "src/lib_alt_QXS/inline_General/*.h")
foreach(header IN LISTS inline_headers)
  get_filename_component(name "${header}" NAME)
  if(NOT name STREQUAL "define_vlen.h")
    configure_file("${header}" "${qxs_inline}/${name}" COPYONLY)
  endif()
endforeach()
configure_file(src/lib_alt_QXS/inline_General/define_vlen.h "${generated}/upstream_vlen.h" COPYONLY)
file(READ "${generated}/upstream_vlen.h" vlen)
# 24^3 needs FP32 2x2 packing with 4x4 SAP blocks in XY.
string(REGEX REPLACE "(#define VLENS[ \t]+)16" "\\14" vlen "${vlen}")
string(REGEX REPLACE "(#define VLEN[XY]S[ \t]+)4" "\\12" vlen "${vlen}")
file(CONFIGURE OUTPUT "${qxs_inline}/define_vlen.h" CONTENT "${vlen}" @ONLY)

file(GLOB_RECURSE bridge_sources CONFIGURE_DEPENDS
  "src/lib/*.cpp" "src/lib_alt/*.cpp" "src/lib_alt_QXS/*.cpp")
list(FILTER bridge_sources EXCLUDE REGEX
  "/(Communicator/(Single|FJMPI)|ResourceManager/OpenMP_stub|(Field|Force/Fermion|Measurements/Fermion)/Org|inline_ACLE)/")
add_library(bridge STATIC ${bridge_sources})
target_include_directories(bridge PUBLIC "${generated}"
  src src/lib src/lib_alt src/lib_alt_QXS)
target_compile_definitions(bridge PUBLIC PC_GNU USE_MPI USE_IMP USE_GROUP_SU3
  USE_OPENMP USE_STD_COMPLEX LIB_CPP11 USE_FACTORY USE_ALT_CODE USE_ALT_QXS
  OMPI_SKIP_MPICXX)
target_link_libraries(bridge PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX m)

add_executable(clover-solve solver/clover_solve.cpp)
target_link_libraries(clover-solve PRIVATE bridge)

add_executable(bridge-hmc sample_hmc/sample_HMC_Wilson_Leapfrog_Nf2.cpp)
target_link_libraries(bridge-hmc PRIVATE bridge)
