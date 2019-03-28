# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( cgal-nef )
set(CXX_STANDARD 14)


cmake_minimum_required(VERSION 2.8.11)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core)

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package

set(CXX_STANDARD 14)

# Creating entries for target: cgal-nef
# ############################

add_executable( cgal-nef  cgal-nef.cpp )

add_to_cached_list( CGAL_EXECUTABLE_TARGETS cgal-nef )

# Link the executable to CGAL and third-party libraries
target_link_libraries(cgal-nef   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
target_compile_features(cgal-nef PUBLIC cxx_std_14)

