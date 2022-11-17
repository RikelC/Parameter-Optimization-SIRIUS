#------------------------------------ 
# CMake file for Parameter Optimization program
#
# Author: Rikel CHAKMA
# Date: 24/06/2022
#-------------------------------------
cmake_minimum_required(VERSION 2.8.12)
project(optimizeParameterProject)
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)
#---------------------
#  C++ compiler
#---------------------
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++ support. Use another C++ compiler.")
endif()

set(CMAKE_C_FLAGS "-O3 -g")
set(CMAKE_CXX_FLAGS "-O3 -g -Wall -Wno-overloaded-virtual -Wno-unknown-pragmas -Wno-write-strings -Wno-deprecated -fpermissive -fPIC -std=c++11")
#----------------------------------------------------------
# Cmake install prefix = install directory used by install
#-----------------------------------------------------------
if(CMAKE_INSTALL_PREFIX)
   message (STATUS "CMAKE_INSTALL_PREFIX found : ${CMAKE_INSTALL_PREFIX} ")
else()
   message (WARNING "CMAKE_INSTALL_PREFIX not found !!! ")
endif()
#-----------------------
# User Library
#--------------------------
set(UserLib_Path /home/sirius/Chakma/Gru_dev/GRU_SIRIUS/GUser-sources/UserLib/build)
find_library(User_LIBRARIES NAMES UserLIB  HINTS ${UserLib_Path}) 

if (User_LIBRARIES)
	message(STATUS "-- found  ${User_LIBRARIES}")
	else()
	message(FATAL_ERROR "-- User_LIBRARIES not found in  ${User_LIBRARIES}")
endif()

set (UserLibIncludes ${UserLib_Path}/../includes)
set (UserLibSources ${UserLib_Path}/../sources)

include_directories(${UserLibIncludes} ${UserLibSources})

#add_subdirectory(UserLib)


#----------------------------------
# ROOT
#---------------------------------
#set(ROOT_DIR /home/global/root/root_v5.34.38_C7)
#list(APPEND CMAKE_PREFIX_PATH ${ROOT_DIR})
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
if(ROOT_FOUND)
	message(STATUS "ROOT found ...")
else()
	message(FATAL_ERROR "!!! ROOT not found. !!! Aborting...")
endif()

include_directories(${ROOT_INCLUDE_DIRS})
#set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS}")
#set(CMAKE_C_FLAGS "${ROOT_CC_FLAGS}")

#----------------------------
# Source and header files
#--------------------------
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/sources/*.cc)

file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/includes/*.hh)
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/includes)

#--------------------------
# Create shared object
#-------------------------
add_executable(optimizeParameter main.cc ${SOURCES} ${HEADERS})
target_link_libraries(optimizeParameter ${ROOT_LIBRARIES} ${User_LIBRARIES} Spectrum)

add_custom_target(optimizeParameterProject DEPENDS optimizeParameter)


install(TARGETS optimizeParameter  DESTINATION bin)
