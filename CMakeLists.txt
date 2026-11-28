cmake_minimum_required(VERSION 3.18)
project(Bridge213 VERSION 2.1.3 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(OpenMP REQUIRED COMPONENTS CXX)
include(CTest)

set(CLOVER_QXS "general" CACHE STRING "CPU kernel backend: general, avx512 or sve512")
set_property(CACHE CLOVER_QXS PROPERTY STRINGS general avx512 sve512)
set(kernel_headers portable)
if(CLOVER_QXS STREQUAL "sve512")
  set(kernel_headers sve512)
  set(kernel_flags -march=armv8.2-a+sve+sve2 -msve-vector-bits=512)
  set(kernel_probe [=[
#include <arm_sve.h>
static_assert(__ARM_FEATURE_SVE_BITS == 512, "Expected 512-bit SVE");
int main() { return 0; }
]=])
elseif(CLOVER_QXS STREQUAL "avx512")
  set(kernel_flags -mavx512f -mprefer-vector-width=512)
  set(kernel_probe [=[
#include <immintrin.h>
#ifndef __AVX512F__
#error AVX-512F required
#endif
int main() { return 0; }
]=])
elseif(NOT CLOVER_QXS STREQUAL "general")
  message(FATAL_ERROR "Use -DCLOVER_QXS=general, avx512 or sve512")
endif()
if(kernel_flags)
  include(CheckCXXSourceCompiles)
  string(JOIN " " CMAKE_REQUIRED_FLAGS ${kernel_flags})
  check_cxx_source_compiles("${kernel_probe}" CLOVER_HAS_${CLOVER_QXS}_COMPILER)
  unset(CMAKE_REQUIRED_FLAGS)
  if(NOT CLOVER_HAS_${CLOVER_QXS}_COMPILER)
    message(FATAL_ERROR "Compiler does not support ${CLOVER_QXS} kernel flags: ${kernel_flags}")
  endif()
endif()

set(generated "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(CONFIGURE OUTPUT "${generated}/version.h"
  CONTENT "#define BRIDGE_VERSION \"${PROJECT_VERSION}\"\n")

file(GLOB_RECURSE bridge_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB kernel_sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/kernels/*.cpp")
if(kernel_flags)
  set_source_files_properties(${kernel_sources} PROPERTIES COMPILE_OPTIONS "${kernel_flags}")
endif()
set(qxs_precision_sources)
foreach(name afield afopr_Clover asolver_SAP_MINRES bridgeQXS_Clover bridgeQXS_Wilson)
  list(APPEND qxs_precision_sources "${CMAKE_CURRENT_SOURCE_DIR}/src/kernels/${name}.cpp")
endforeach()
list(REMOVE_ITEM bridge_sources ${qxs_precision_sources})

add_library(bridge_config INTERFACE)
target_include_directories(bridge_config INTERFACE "${generated}" src "src/kernels/${kernel_headers}")
target_compile_definitions(bridge_config INTERFACE PC_GNU USE_MPI USE_GROUP_SU3
  USE_OPENMP USE_STD_COMPLEX LIB_CPP11 USE_FACTORY USE_ALT_CODE USE_ALT_QXS
  OMPI_SKIP_MPICXX)
target_link_libraries(bridge_config INTERFACE MPI::MPI_CXX OpenMP::OpenMP_CXX m)

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
