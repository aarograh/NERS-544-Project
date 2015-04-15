# AUTHORS: Aaron Graham, Mike Jarrett
# PURPOSE: NERS 544 Course Project
# DATE   : April 3, 2015

# Basic setup stuff
cmake_minimum_required(VERSION 2.8)

# Set path for cmake
set(CMAKE_MODULE_PATH
  "${CMAKE_SOURCE_DIR}/cmake"
  ${CMAKE_MODULE_PATH}
)

# Define project settings
project(NERS-544-Project)
set(APPLICATION_NAME ${PROJECT_NAME})
IF(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Includes
include(CMakeCompilerFlags)
include(MacroEnsureOutOfSourceBuild)
include(CMakeDefaults)

# Check paths
MACRO_ENSURE_OUT_OF_SOURCE_BUILD("Attempted to build in source directory!")

# Compile Program and Docs
include_directories(inc)
add_subdirectory(src)
add_subdirectory(doc)
add_executable(NERS-544.exe ${SOURCE} ${HEADERS})
