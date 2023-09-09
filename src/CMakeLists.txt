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

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, PROFILING, and TESTING.  
INCLUDE(SetFortranFlags) 

# Add the source files
SET(SRC_FILES
    ${SRC}/globals/module_globals.f90
    ${SRC}/util/module_util.f90
    ${SRC}/io/module_io.f90
    ${SRC}/crater/module_crater.f90
    ${SRC}/ejecta/module_ejecta.f90
    ${SRC}/init/module_init.f90
    ${SRC}/seismic/module_seismic.f90
    ${SRC}/regolith/module_regolith.f90
    ${SRC}/porosity/module_porosity.f90
    ${SRC}/realistic/module_realistic.f90
    ${SRC}/util/util_toupper.f90
    ${SRC}/util/util_mrgrnk.f90
    ${SRC}/util/util_add_to_layer.f90
    ${SRC}/util/util_remove_from_layer.f90
    ${SRC}/util/util_sort_layer.f90
    ${SRC}/util/util_periodic.f90
    ${SRC}/util/util_search.f90
    ${SRC}/util/util_rootbracketed.f90
    ${SRC}/util/util_diffusion_solver.f90
    ${SRC}/util/util_area_intersection.f90
    ${SRC}/util/util_poisson.f90
    ${SRC}/util/util_pop.f90
    ${SRC}/util/util_pop_array.f90
    ${SRC}/util/util_push_array.f90
    ${SRC}/util/util_traverse_pop_array.f90
    ${SRC}/util/util_destroy_list.f90
    ${SRC}/util/util_init_array.f90
    ${SRC}/util/util_perlin_noise.f90
    ${SRC}/util/util_random_number_normal.f90
    ${SRC}/util/util_npf_timefunc.f90
    ${SRC}/util/util_tscale.f90
    ${SRC}/util/util_t_from_scale.f90
    ${SRC}/io/io_read_const.f90
    ${SRC}/io/io_get_token.f90
    ${SRC}/io/io_input.f90
    ${SRC}/io/io_read_craterlist.f90
    ${SRC}/io/io_read_prod.f90
    ${SRC}/io/io_read_regotrack.f90
    ${SRC}/io/io_read_vdist.f90
    ${SRC}/io/io_read_surf.f90
    ${SRC}/io/io_ejecta_table.f90
    ${SRC}/io/io_write_dist.f90
    ${SRC}/io/io_write_tally.f90
    ${SRC}/io/io_write_surf.f90
    ${SRC}/io/io_write_const.f90
    ${SRC}/io/io_write_regotrack.f90
    ${SRC}/io/io_crater_profile.f90
    ${SRC}/io/io_updatePbar.f90
    ${SRC}/io/io_resetPbar.f90
    ${SRC}/io/io_splash.f90
    ${SRC}/io/io_write_pindex_map.f90
    ${SRC}/ejecta/ejecta_emplace.f90
    ${SRC}/ejecta/ejecta_ray_pattern.f90
    ${SRC}/ejecta/ejecta_ray_pattern_func.f90
    ${SRC}/ejecta/ejecta_ray_func.f90
    ${SRC}/ejecta/ejecta_blanket.f90
    ${SRC}/ejecta/ejecta_blanket_func.f90
    ${SRC}/ejecta/ejecta_thickness.f90
    ${SRC}/ejecta/ejecta_table_define.f90
    ${SRC}/ejecta/ejecta_interpolate.f90
    ${SRC}/ejecta/ejecta_rootfind.f90
    ${SRC}/ejecta/ejecta_soften.f90
    ${SRC}/ejecta/ejecta_distance_estimate.f90
    ${SRC}/crater/crater_populate.f90
    ${SRC}/crater/crater_scale.f90
    ${SRC}/crater/crater_generate.f90
    ${SRC}/crater/crater_dimensions.f90
    ${SRC}/crater/crater_averages.f90
    ${SRC}/crater/crater_emplace.f90
    ${SRC}/crater/crater_realistic_topography.f90
    ${SRC}/crater/crater_form_interior.f90
    ${SRC}/crater/crater_form_exterior.f90
    ${SRC}/crater/crater_form_exterior_func.f90
    ${SRC}/crater/crater_form_exterior_rootfind.f90
    ${SRC}/crater/crater_record.f90
    ${SRC}/crater/crater_tally_true.f90
    ${SRC}/crater/crater_tally_observed.f90
    ${SRC}/crater/crater_slope_collapse.f90
    ${SRC}/crater/crater_soften.f90
    ${SRC}/crater/crater_subpixel_diffusion.f90
    ${SRC}/crater/crater_make_list.f90
    ${SRC}/crater/crater_degradation_function.f90
    ${SRC}/crater/crater_visibility.f90
    ${SRC}/crater/crater_get_degradation_state.f90
    ${SRC}/crater/crater_profile.f90
    ${SRC}/crater/crater_superdomain.f90
    ${SRC}/init/init_domain.f90
    ${SRC}/init/init_dist.f90
    ${SRC}/init/init_surf.f90
    ${SRC}/init/init_regolith_stack.f90
    ${SRC}/seismic/seismic_shake.f90
    ${SRC}/seismic/seismic_distance.f90
    ${SRC}/seismic/seismic_kdiff_func.f90
    ${SRC}/regolith/regolith_circle_sector_func.f90
    ${SRC}/regolith/regolith_cubic_func.f90
    ${SRC}/regolith/regolith_melt_func.f90
    ${SRC}/regolith/regolith_melt_fraction.f90
    ${SRC}/regolith/regolith_melt_zone.f90
    ${SRC}/regolith/regolith_quadratic_func.f90
    ${SRC}/regolith/regolith_quartic_func.f90
    ${SRC}/regolith/regolith_mix.f90
    ${SRC}/regolith/regolith_depth_model.f90
    ${SRC}/regolith/regolith_streamtube.f90
    ${SRC}/regolith/regolith_streamtube_head.f90
    ${SRC}/regolith/regolith_streamtube_lineseg.f90
    ${SRC}/regolith/regolith_subpixel_streamtube.f90
    ${SRC}/regolith/regolith_transport.f90
    ${SRC}/regolith/regolith_traverse_streamtube.f90
    ${SRC}/regolith/regolith_subcrater_mix.f90
    ${SRC}/regolith/regolith_melt_glass.f90
    ${SRC}/regolith/regolith_superdomain.f90
    ${SRC}/regolith/regolith_melt_zone_superdomain.f90
    ${SRC}/regolith/regolith_streamtube_volume_func.f90
    ${SRC}/regolith/regolith_shock_damage_zone.f90
    ${SRC}/regolith/regolith_shock_damage.f90
    ${SRC}/regolith/regolith_interior.f90
    ${SRC}/porosity/porosity_form_interior.f90
    ${SRC}/realistic/realistic_perlin_noise.f90
    ${SRC}/realistic/realistic_crater_topography.f90
)


SET(DRIVER_src ${SRC}/main/CTEM.f90)

# Combine all source files 
set(CTEM_src ${SRC_FILES})

# Turn preprocessor on for all files
SET_SOURCE_FILES_PROPERTIES(${CTEM_src} ${DRIVER_src} PROPERTIES Fortran_PREPROCESS ON)

# Set Build Type variable
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "DEBUG")
   ADD_DEFINITIONS(-DDEBUG)
ELSEIF(BT STREQUAL "PROFILE")
   ADD_DEFINITIONS(-DPROFILE)
ENDIF()

#####################################################
# Add the needed libraries 
#####################################################
# Create a library from the source files, except the driver
ADD_LIBRARY(${CTEM_LIBRARY} ${CTEM_src})
IF (NOT BUILD_SHARED_LIBS) 
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} PROPERTY POSITION_INDEPENDENT_CODE)
ENDIF ()


# Define the executable name
SET(CTEM_DRIVER ctem_exe)
ADD_EXECUTABLE(${CTEM_DRIVER} ${DRIVER_src})
SET_PROPERTY(TARGET ${CTEM_DRIVER} PROPERTY OUTPUT_NAME ctem)

IF(USE_OPENMP OR USE_SIMD)
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY COMPILE_FLAGS  "${OpenMP_Fortran_FLAGS} ")
    SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY LINK_FLAGS  "${OpenMP_Fortran_FLAGS} ")
ENDIF()

IF (CMAKE_SYSTEM_NAME STREQUAL "Windows")
   SET_PROPERTY(TARGET ${CTEM_LIBRARY} ${CTEM_DRIVER} APPEND_STRING PROPERTY LINK_FLAGS  "/NODEFAULTLIB")
ENDIF()


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

