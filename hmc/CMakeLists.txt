# This module is also installed by the top-level Clover-SVE CMake project.
if(NOT TARGET bridge)
  message(FATAL_ERROR "Enable the Bridge++ CPU library before adding hmc/")
endif()
set(hmc_runtime "${CMAKE_BINARY_DIR}/libexec/clover-hmc")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/bin" "${hmc_runtime}")
configure_file(clover-hmc "${CMAKE_BINARY_DIR}/bin/clover-hmc" COPYONLY)
foreach(script generator.py engine.py gauge_io.py qdp_stream.py profile.json)
  configure_file("runtime/${script}" "${hmc_runtime}/${script}" COPYONLY)
endforeach()

foreach(kind gauge light logdet strange)
  add_executable(clover-hmc-${kind} native/${kind}.cpp)
  target_link_libraries(clover-hmc-${kind} PRIVATE bridge)
  set_target_properties(clover-hmc-${kind} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${hmc_runtime}")
  target_compile_features(clover-hmc-${kind} PRIVATE cxx_std_11)
  install(TARGETS clover-hmc-${kind} RUNTIME DESTINATION libexec/clover-hmc)
endforeach()
install(PROGRAMS clover-hmc DESTINATION bin)
install(FILES runtime/generator.py runtime/engine.py runtime/gauge_io.py runtime/qdp_stream.py runtime/profile.json DESTINATION libexec/clover-hmc)


if(BUILD_TESTING AND NOT CMAKE_CROSSCOMPILING)
  find_package(Python3 QUIET COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME hmc_gauge_validation COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/gauge_validation.py)
    set_tests_properties(hmc_gauge_validation PROPERTIES SKIP_RETURN_CODE 77)
    add_test(NAME native_chroma_energy_order COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/tests/chroma_energy_order.py ${CMAKE_CURRENT_SOURCE_DIR}/runtime
      ${PROJECT_SOURCE_DIR}/validation/L24-hmc-cfg3000/chroma-energy-order.json)
    set_tests_properties(native_chroma_energy_order PROPERTIES SKIP_RETURN_CODE 77)
  endif()
endif()
