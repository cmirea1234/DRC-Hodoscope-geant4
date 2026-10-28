cmake_minimum_required(VERSION 3.16)
project(HodoscopeKNUTwoModels LANGUAGES CXX)
# Match the tested version. Library patches are documented separately.
find_package(Geant4 11.4.1 EXACT REQUIRED)
include(${Geant4_USE_FILE})
foreach(model IN ITEMS original light32)
  set(sources models/${model}/src/main.cc
              models/${model}/src/DetectorConstruction.cc
              models/${model}/src/ActionInitialization.cc
              models/${model}/src/V21GeometryTolerance.cc)
  if(model STREQUAL "light32")
    list(APPEND sources models/light32/src/StraightFibres.cc)
  endif()
  add_executable(hodoscope-${model} ${sources})
  target_include_directories(hodoscope-${model} PRIVATE models/${model}/include)
  target_compile_features(hodoscope-${model} PRIVATE cxx_std_17)
  target_link_libraries(hodoscope-${model} PRIVATE ${Geant4_LIBRARIES})
endforeach()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
enable_testing()
add_test(NAME portable_bundle_integrity
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/verify_bundle.py)
