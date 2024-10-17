cmake_minimum_required(VERSION 3.16.0)

project(molbar C Fortran)

option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
endif()

if (NOT SKBUILD)
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/../cmake/")
endif()

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

find_package(PythonExtensions REQUIRED)
find_package(NumPy)
find_package(F2PY REQUIRED)

message(STATUS "Cmake build directory: ${CMAKE_BINARY_DIR}$")
message(STATUS "F2PY included from: ${F2PY_INCLUDE_DIR}")

# Osipov-Pickup-Dunmir Index 
set(f2py_module_name "opd")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/molbar/indices/opd.f90")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )

add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${fortran_src_file})

target_include_directories(${generated_module_file} PUBLIC
                           ${F2PY_INCLUDE_DIRS}
                           ${PYTHON_INCLUDE_DIRS})
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")


if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()


if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib)
else()
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib)
endif()

# Bond Order

# Build libary
add_library(unsaturated_graph STATIC "${CMAKE_SOURCE_DIR}/molbar/molgraph/edges/unsaturated_graph.f90")
set_property(TARGET unsaturated_graph PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(unsaturated_edges STATIC "${CMAKE_SOURCE_DIR}/molbar/molgraph/edges/unsaturated_edges.f90")
set_property(TARGET unsaturated_edges PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(valence_conditions STATIC "${CMAKE_SOURCE_DIR}/molbar/molgraph/edges/valence_conditions.f90")
set_property(TARGET valence_conditions PROPERTY POSITION_INDEPENDENT_CODE ON)


set(f2py_module_name "bo_matrix")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/molbar/molgraph/edges/bo_matrix.f90")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )


add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${fortran_src_file})

add_dependencies(${generated_module_file} unsaturated_graph unsaturated_edges valence_conditions)

target_link_libraries(${generated_module_file} unsaturated_graph)
target_link_libraries(${generated_module_file} unsaturated_edges)
target_link_libraries(${generated_module_file} valence_conditions)
target_include_directories(${generated_module_file} PUBLIC
                           ${F2PY_INCLUDE_DIRS}
                           ${PYTHON_INCLUDE_DIRS})
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")


if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()


if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib)
else()
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib)
endif()


# Force field

# Build libary
add_library(fortran_helper STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/fortran_helper.f90")
set_property(TARGET fortran_helper PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(a_angle_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/a_angle_derivatives.f90")
set_property(TARGET a_angle_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(b_angle_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/b_angle_derivatives.f90")
set_property(TARGET b_angle_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(c_angle_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/c_angle_derivatives.f90")
set_property(TARGET c_angle_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(a_dihedral_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/a_dihedral_derivatives.f90")
set_property(TARGET a_dihedral_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(b_dihedral_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/b_dihedral_derivatives.f90")
set_property(TARGET b_dihedral_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(c_dihedral_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/c_dihedral_derivatives.f90")
set_property(TARGET c_dihedral_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(d_dihedral_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/d_dihedral_derivatives.f90")
set_property(TARGET d_dihedral_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(inv_rij_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/inv_rij_derivatives.f90")
set_property(TARGET inv_rij_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(bond_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/bond_derivatives.f90")
set_property(TARGET bond_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(angle_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/angle_derivatives.f90")
set_property(TARGET angle_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(dihedral_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/dihedral_derivatives.f90")
set_property(TARGET dihedral_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)
# build libary
add_library(repulsion_derivatives STATIC "${CMAKE_SOURCE_DIR}/molbar/force_field/repulsion_derivatives.f90")
set_property(TARGET repulsion_derivatives PROPERTY POSITION_INDEPENDENT_CODE ON)

target_link_libraries(a_angle_derivatives fortran_helper)
target_link_libraries(b_angle_derivatives fortran_helper)
target_link_libraries(c_angle_derivatives fortran_helper b_angle_derivatives)
target_link_libraries(a_dihedral_derivatives fortran_helper)
target_link_libraries(b_dihedral_derivatives fortran_helper)
target_link_libraries(c_dihedral_derivatives fortran_helper b_dihedral_derivatives)
target_link_libraries(d_dihedral_derivatives fortran_helper)
target_link_libraries(inv_rij_derivatives fortran_helper)
target_link_libraries(bond_derivatives fortran_helper)
target_link_libraries(angle_derivatives fortran_helper a_angle_derivatives c_angle_derivatives b_angle_derivatives)
target_link_libraries(dihedral_derivatives fortran_helper a_dihedral_derivatives c_dihedral_derivatives b_dihedral_derivatives  d_dihedral_derivatives)
target_link_libraries(repulsion_derivatives fortran_helper inv_rij_derivatives)

target_link_libraries(angle_derivatives fortran_helper a_angle_derivatives c_angle_derivatives b_angle_derivatives)
target_link_libraries(bond_derivatives fortran_helper b_angle_derivatives)
target_link_libraries(dihedral_derivatives fortran_helper a_dihedral_derivatives c_dihedral_derivatives b_dihedral_derivatives  d_dihedral_derivatives)
target_link_libraries(repulsion_derivatives fortran_helper inv_rij_derivatives)

set(f2py_module_name "analytical_derivatives")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/molbar/force_field/analytical_derivatives.f90")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )

add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${fortran_src_file})

target_link_libraries(${generated_module_file} bond_derivatives)
target_link_libraries(${generated_module_file} angle_derivatives)
target_link_libraries(${generated_module_file} dihedral_derivatives)
target_link_libraries(${generated_module_file} repulsion_derivatives)
target_include_directories(${generated_module_file} PUBLIC
                           ${F2PY_INCLUDE_DIRS}
                           ${PYTHON_INCLUDE_DIRS})
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")


if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()


if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib/)
else()
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib/)
endif()


# Force field energy
# Build libary
set(f2py_module_name "energy")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/molbar/force_field/energy.f90")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${F2PY_EXECUTABLE}
    -m ${f2py_module_name}
    ${fortran_src_file}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )

add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${fortran_src_file})

add_dependencies(${generated_module_file} fortran_helper)

target_link_libraries(${generated_module_file} fortran_helper)

target_include_directories(${generated_module_file} PUBLIC
                           ${F2PY_INCLUDE_DIRS}
                           ${PYTHON_INCLUDE_DIRS})
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")


if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()


if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib/)
else()
  install(TARGETS ${generated_module_file} DESTINATION molbar/fortranlib/)
endif()
