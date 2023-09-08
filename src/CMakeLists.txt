# Copyright 2023 - David Minton, Carlisle Wishard, Jennifer Pouplin, Jake Elliott, & Dana Singh
# This file is part of Swiftest.
# Swiftest is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# Swiftest is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with Swiftest. 
# If not, see: https://www.gnu.org/licenses. 

# Communicate version number and other CMake build variables to the source code
#set(GLOBAL_MODULE_IN ${SRC}/globals/globals_module.f90.in)
#set(GLOBAL_MODULE_OUT ${SRC}/globals/globals_module.f90)
#CONFIGURE_FILE(${GLOBAL_MODULE_IN} ${GLOBAL_MODULE_OUT})

IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

INCLUDE(SetParallelizationLibrary)
INCLUDE(SetUpNetCDF)
IF (COMPILER_OPTIONS STREQUAL "Intel" AND NOT CMAKE_SYSTEM_NAME STREQUAL "Windows")
    INCLUDE(SetMKL)
ENDIF ()

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, PROFILING, and TESTING.  
INCLUDE(SetFortranFlags) 

# Add the source files
SET(SOURCE_FILES

)


SET(DRIVER_src ${SRC}/main/CTEM.f90)

# Combine all source files 
set(CTEM_src ${SOURCE_FILES})

# Turn preprocessor on for all files
SET_SOURCE_FILES_PROPERTIES(${CTEM_src} ${DRIVER_src} PROPERTIES Fortran_PREPROCESS ON)

#Set strict vs fast math flags
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
IF(BT STREQUAL "RELEASE" OR BT STREQUAL "PROFILE")
   SET_PROPERTY(SOURCE ${STRICT_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${STRICTMATH_FLAGS}")
   SET_PROPERTY(SOURCE ${FAST_MATH_FILES} APPEND_STRING PROPERTY COMPILE_FLAGS "${FASTMATH_FLAGS}")
ENDIF()

IF(BT STREQUAL "DEBUG")
   ADD_DEFINITIONS(-DDEBUG)
ELSEIF(BT STREQUAL "PROFILE")
   ADD_DEFINITIONS(-DPROFILE)
ENDIF()

# Define the executable name
SET(CTEM_DRIVER ctem_driver)
ADD_EXECUTABLE(${CTEM_DRIVER} ${DRIVER_src})

#####################################################
# Add the needed libraries 
#####################################################
# Create a library from the source files, except the driver
ADD_LIBRARY(${CTEM_LIBRARY} ${CTEM_src})
IF (NOT BUILD_SHARED_LIBS) 
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} PROPERTY POSITION_INDEPENDENT_CODE)
ENDIF ()

TARGET_INCLUDE_DIRECTORIES(${CTEM_LIBRARY} PUBLIC ${NETCDF_INCLUDE_DIR})
TARGET_INCLUDE_DIRECTORIES(${CTEM_DRIVER} PUBLIC ${NETCDF_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(${CTEM_LIBRARY} PUBLIC ${NETCDF_LIBRARIES})
TARGET_LINK_LIBRARIES(${CTEM_DRIVER} PUBLIC ${CTEM_LIBRARY} ${NETCDF_LIBRARIES})

IF(USE_OPENMP OR USE_SIMD)
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY LINK_FLAGS  "${OpenMP_Fortran_FLAGS} ")
ENDIF()

IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
   SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY LINK_FLAGS  "/NODEFAULTLIB")
ENDIF()


IF(USE_COARRAY)
    TARGET_COMPILE_DEFINITIONS(${CTEM_LIBRARY} PUBLIC -DCOARRAY)
    TARGET_COMPILE_DEFINITIONS(${CTEM_DRIVER} PUBLIC -DCOARRAY)
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY COMPILE_FLAGS  "${Coarray_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY LINK_FLAGS  "${Coarray_Fortran_FLAGS} ")
ENDIF(USE_COARRAY)

# Check to see if the compiler allows for local-spec in do concurrent statements. Set a preprocessor variable if it does
SET(TESTFILE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}")
SET(TESTFILE "${TESTFILE}/CMakeTmp/testFortranDoConcurrentLoc.f90")
FILE(WRITE "${TESTFILE}"
"
program TestDoConcurrentLoc
integer :: i
real,dimension(10) :: a
do concurrent(i = 1:10) shared(a)
   a(i) = i 
end do
end program TestDoConcurrentLoc
")
TRY_COMPILE(DOCONLOC_WORKS ${CMAKE_BINARY_DIR} ${TESTFILE} COMPILE_DEFINITIONS "${CMAKE_Fortran_FLAGS}" OUTPUT_VARIABLE OUTPUT)
IF (DOCONLOC_WORKS)
    MESSAGE(STATUS "DO CONCURRENT supports locality-spec")
    TARGET_COMPILE_DEFINITIONS(${CTEM_DRIVER} PRIVATE -DDOCONLOC)
ELSE ()
    MESSAGE(STATUS "DO CONCURRENT does not support locality-spec")
ENDIF (DOCONLOC_WORKS)

# Check to see if quad precision is supported
SET(TESTFILE "${CMAKE_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}")
SET(TESTFILE "${TESTFILE}/CMakeTmp/testFortranQuadPrecisionReal.f90")
FILE(WRITE "${TESTFILE}"
"
program TestQuadPrecisionReal
integer, parameter :: QP = selected_Real_kind(30) 
real(QP) :: x
end program TestQuadPrecisionReal
")
TRY_COMPILE(QUADPREC ${CMAKE_BINARY_DIR} ${TESTFILE} COMPILE_DEFINITIONS "${CMAKE_Fortran_FLAGS}" OUTPUT_VARIABLE OUTPUT)
IF (QUADPREC)
    MESSAGE(STATUS "Quad precision real is supported")
    TARGET_COMPILE_DEFINITIONS(${CTEM_DRIVER} PRIVATE -DQUADPREC)
ELSE ()
    MESSAGE(STATUS "Quad precision real is not supported")
ENDIF ()


#####################################
# Tell how to install this executable
#####################################
IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
    FILE(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX} CMAKE_INSTALL_PREFIX)
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)
INSTALL(TARGETS ${CTEM_DRIVER} ${CTEM_LIBRARY}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    INCLUDES DESTINATION include
)

