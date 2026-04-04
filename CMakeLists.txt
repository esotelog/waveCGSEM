##
#  CMake script for the step-3 tutorial program:
##

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#    file(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#    file(GLOB_RECURSE TARGET_INC  "include/*.h")
#    set(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.

# Usually, you will not need to modify anything beyond this point...

cmake_minimum_required(VERSION 3.13.4)

find_package(deal.II 9.7.0
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()

# setting up cmake w. deal.II
DEAL_II_INITIALIZE_CACHED_VARIABLES()


#project folder: CMakeLists.txt is in this folder
PROJECT(waveCGSEM LANGUAGES CXX)

#Include and source directories
SET(WAVECGSEM_SOURCES
    source/main_waveCGSEM.cc
    source/functions.cc
    source/Param_handler.cc
    # ADD ALL your project's implementation (.cc) files here!
)

ADD_EXECUTABLE(waveCGSEM_exe ${WAVECGSEM_SOURCES})

TARGET_INCLUDE_DIRECTORIES(waveCGSEM_exe PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

TARGET_INCLUDE_DIRECTORIES(waveCGSEM_exe PUBLIC 
    ${DEAL_II_INCLUDE_DIRS}
)

#Configures exec with deal.II library
DEAL_II_SETUP_TARGET(waveCGSEM_exe)

SET_PROPERTY(TARGET waveCGSEM_exe PROPERTY CXX_STANDARD 17)
SET_PROPERTY(TARGET waveCGSEM_exe PROPERTY CXX_STANDARD_REQUIRED ON)

