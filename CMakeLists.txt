cmake_minimum_required(VERSION 3.18)
project(Bridge213 VERSION 2.1.3 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED COMPONENTS CXX)
include(CTest)

set(CLOVER_QXS "general" CACHE STRING "QXS backend: general or sve512")
set_property(CACHE CLOVER_QXS PROPERTY STRINGS general sve512)
set(qxs_inline_source inline_General)
if(CLOVER_QXS STREQUAL "sve512")
  include(CheckCXXSourceCompiles)
  set(qxs_sve_flags -march=armv8.2-a+sve+sve2 -msve-vector-bits=512)
  string(JOIN " " CMAKE_REQUIRED_FLAGS ${qxs_sve_flags})
  check_cxx_source_compiles([=[
#include <arm_sve.h>
static_assert(__ARM_FEATURE_SVE_BITS == 512, "Expected 512-bit SVE");
int main() { return 0; }
]=] CLOVER_HAS_SVE512_COMPILER)
  unset(CMAKE_REQUIRED_FLAGS)
  if(NOT CLOVER_HAS_SVE512_COMPILER)
    message(FATAL_ERROR "sve512 requires an AArch64 compiler with 512-bit SVE2 support")
  endif()
  set(qxs_inline_source inline_ACLE)
elseif(NOT CLOVER_QXS STREQUAL "general")
  message(FATAL_ERROR "Use -DCLOVER_QXS=general or sve512")
endif()

# Generate the upstream configuration headers inside the build tree.
set(generated "${CMAKE_CURRENT_BINARY_DIR}/generated")
set(bridge_version "${PROJECT_VERSION}")
file(CONFIGURE OUTPUT "${generated}/version.h"
  CONTENT "#define BRIDGE_VERSION \"${bridge_version}\"\n")

set(qxs_inline "${generated}/lib_alt_QXS/inline")
file(GLOB inline_headers CONFIGURE_DEPENDS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "src/lib_alt_QXS/${qxs_inline_source}/*.h" "src/lib_alt_QXS/inline_common/*.h")
foreach(header IN LISTS inline_headers)
  get_filename_component(name "${header}" NAME)
  if(NOT name STREQUAL "define_vlen.h")
    configure_file("${header}" "${qxs_inline}/${name}" COPYONLY)
  endif()
endforeach()
configure_file(src/lib_alt_QXS/${qxs_inline_source}/define_vlen.h "${generated}/upstream_vlen.h" COPYONLY)
file(READ "${generated}/upstream_vlen.h" vlen)
if(CLOVER_QXS STREQUAL "general")
  # 24^3 needs FP32 2x2 packing with 4x4 SAP blocks in XY.
  string(REGEX REPLACE "(#define VLENS[ \t]+)16" "\\14" vlen "${vlen}")
  string(REGEX REPLACE "(#define VLEN[XY]S[ \t]+)4" "\\12" vlen "${vlen}")
endif()
file(CONFIGURE OUTPUT "${qxs_inline}/define_vlen.h" CONTENT "${vlen}" @ONLY)

file(GLOB_RECURSE bridge_sources CONFIGURE_DEPENDS
  "src/lib/*.cpp" "src/lib_alt/*.cpp" "src/lib_alt_QXS/*.cpp")
set(qxs_precision_sources
  src/lib_alt_QXS/BridgeQXS/bridgeQXS_Clover.cpp
  src/lib_alt_QXS/BridgeQXS/bridgeQXS_Wilson.cpp
  src/lib_alt_QXS/Field/afield.cpp
  src/lib_alt_QXS/Fopr/afopr_Clover.cpp
  src/lib_alt_QXS/Solver/asolver_SAP_MINRES.cpp)
foreach(source IN LISTS qxs_precision_sources)
  list(REMOVE_ITEM bridge_sources "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
endforeach()

add_library(bridge_config INTERFACE)
target_include_directories(bridge_config INTERFACE "${generated}" src src/lib src/lib_alt src/lib_alt_QXS)
target_compile_definitions(bridge_config INTERFACE PC_GNU USE_MPI USE_GROUP_SU3
  USE_OPENMP USE_STD_COMPLEX LIB_CPP11 USE_FACTORY USE_ALT_CODE USE_ALT_QXS
  OMPI_SKIP_MPICXX)
target_link_libraries(bridge_config INTERFACE MPI::MPI_CXX OpenMP::OpenMP_CXX m)
if(CLOVER_QXS STREQUAL "sve512")
  target_compile_options(bridge_config INTERFACE ${qxs_sve_flags})
endif()

add_library(bridge STATIC ${bridge_sources})
target_link_libraries(bridge PUBLIC bridge_config)
foreach(precision IN ITEMS 0 1)
  add_library(qxs_precision_${precision} OBJECT ${qxs_precision_sources})
  target_link_libraries(qxs_precision_${precision} PRIVATE bridge_config)
  target_compile_definitions(qxs_precision_${precision} PRIVATE QXS_SINGLE_PRECISION=${precision})
  target_sources(bridge PRIVATE $<TARGET_OBJECTS:qxs_precision_${precision}>)
endforeach()

add_executable(clover-solve solver/clover_solve.cpp)
target_link_libraries(clover-solve PRIVATE bridge)

add_subdirectory(hmc)
install(TARGETS clover-solve RUNTIME DESTINATION bin)

if(BUILD_TESTING AND NOT CMAKE_CROSSCOMPILING)
  add_executable(clover-base-check tests/base_operator.cpp)
  target_link_libraries(clover-base-check PRIVATE bridge)
  add_executable(clover-sap-check tests/block_check.cpp)
  target_link_libraries(clover-sap-check PRIVATE bridge)
  add_test(NAME analytic_clover_plane_wave COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:clover-base-check> ${MPIEXEC_POSTFLAGS})
  add_test(NAME sap_half_pack_operator COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:clover-sap-check> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(analytic_clover_plane_wave sap_half_pack_operator PROPERTIES
    ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 120)
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME numa_rank_bindings COMMAND ${Python3_EXECUTABLE}
      "${PROJECT_SOURCE_DIR}/tests/numa_rank_bindings.py")
    foreach(test solver_layouts solver_checkpoint solver_gcr)
      add_test(NAME ${test} COMMAND ${Python3_EXECUTABLE} "${PROJECT_SOURCE_DIR}/tests/${test}.py"
        ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:clover-solve> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(${test} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" TIMEOUT 180)
    endforeach()
  endif()
endif()
