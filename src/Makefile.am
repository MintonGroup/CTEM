bin_PROGRAMS = CTEM 
#ifort optimized flags
SIMDVEC = -simd -xhost -align all -assume contiguous_assumed_shape -vecabi=cmdtarget -fp-model no-except -fma
PAR = -qopenmp -parallel  
HEAPARR = -heap-arrays 
OPTREPORT = -qopt-report=5 
IPRODUCTION = -g -traceback -no-wrap-margin -assume byterecl -O3 -qopt-prefetch=0 -sox $(PAR) $(SIMDVEC) $(HEAPARR)
#IDEBUG = -O0 -g -traceback -debug all -nogen-interfaces -assume byterecl -m64 -heap-arrays -FR -no-pie -no-ftz -fpe-all=0 -mp1 -fp-model strict -fpe0 -align all -pad -ip -prec-div -prec-sqrt -assume protect-parens -CB -no-wrap-margin -init=snan,arrays 
IDEBUG = -O0 -g -traceback -debug all -nogen-interfaces -assume byterecl -m64 -heap-arrays -FR -no-pie -no-ftz -fpe-all=0 -mp1 -fp-model strict -fpe0 -align all -pad -ip -prec-sqrt -assume protect-parens -CB -no-wrap-margin -init=snan,arrays 
#AM_FCFLAGS =  $(IPRODUCTION)
#ifort debug flags

#gfortran optimized flags
AM_FCFLAGS =  -O3 -fopenmp -ffree-form  -g -fbounds-check -fbacktrace
#gfortran debug flags
#AM_FCFLAGS =  -O0 -g -fopenmp -fbounds-check -Wall -Warray-bounds -Warray-temporaries -Wimplicit-interface -ffree-form -fsanitize-address-use-after-scope -fstack-check  -fsanitize=bounds-strict -fsanitize=undefined -fsanitize=signed-integer-overflow -fsanitize=object-size -fstack-protector-all

CTEM_SOURCES = globals/module_globals.f90\
util/module_util.f90\
io/module_io.f90\
crater/module_crater.f90\
ejecta/module_ejecta.f90\
init/module_init.f90\
seismic/module_seismic.f90\
regolith/module_regolith.f90\
porosity/module_porosity.f90\
realistic/module_realistic.f90\
util/util_toupper.f90\
util/util_mrgrnk.f90\
util/util_add_to_layer.f90\
util/util_remove_from_layer.f90\
util/util_sort_layer.f90\
util/util_periodic.f90\
util/util_search.f90\
util/util_rootbracketed.f90\
util/util_diffusion_solver.f90\
util/util_area_intersection.f90\
util/util_poisson.f90\
util/util_pop_array.f90\
util/util_push_array.f90\
util/util_traverse_pop_array.f90\
util/util_destroy_list.f90\
util/util_init_array.f90\
util/util_perlin_noise.f90\
io/io_read_const.f90\
io/io_get_token.f90\
io/io_input.f90\
io/io_read_craterlist.f90\
io/io_read_prod.f90\
io/io_read_regotrack.f90\
io/io_read_vdist.f90\
io/io_read_surf.f90\
io/io_ejecta_table.f90\
io/io_write_dist.f90\
io/io_write_tally.f90\
io/io_write_surf.f90\
io/io_write_const.f90\
io/io_write_regotrack.f90\
io/io_crater_profile.f90\
io/io_updatePbar.f90\
io/io_resetPbar.f90\
io/io_splash.f90\
io/io_write_pindex_map.f90\
ejecta/ejecta_emplace.f90\
ejecta/ejecta_ray_pattern.f90\
ejecta/ejecta_ray_pattern_func.f90\
ejecta/ejecta_ray_func.f90\
ejecta/ejecta_blanket.f90\
ejecta/ejecta_blanket_func.f90\
ejecta/ejecta_thickness.f90\
ejecta/ejecta_table_define.f90\
ejecta/ejecta_interpolate.f90\
ejecta/ejecta_rootfind.f90\
ejecta/ejecta_soften.f90\
ejecta/ejecta_distance_estimate.f90\
crater/crater_populate.f90\
crater/crater_scale.f90\
crater/crater_generate.f90\
crater/crater_dimensions.f90\
crater/crater_averages.f90\
crater/crater_emplace.f90\
crater/crater_realistic_topography.f90\
crater/crater_form_interior.f90\
crater/crater_form_exterior.f90\
crater/crater_form_exterior_func.f90\
crater/crater_form_exterior_rootfind.f90\
crater/crater_record.f90\
crater/crater_tally_true.f90\
crater/crater_tally_observed.f90\
crater/crater_slope_collapse.f90\
crater/crater_soften.f90\
crater/crater_subpixel_diffusion.f90\
crater/crater_make_list.f90\
crater/crater_degradation_function.f90\
crater/crater_visibility.f90\
crater/crater_get_degradation_state.f90\
crater/crater_profile.f90\
crater/crater_superdomain.f90\
init/init_domain.f90\
init/init_dist.f90\
init/init_surf.f90\
init/init_regolith_stack.f90\
seismic/seismic_shake.f90\
seismic/seismic_distance.f90\
seismic/seismic_kdiff_func.f90\
regolith/regolith_circle_sector_func.f90\
regolith/regolith_cubic_func.f90\
regolith/regolith_melt_func.f90\
regolith/regolith_melt_fraction.f90\
regolith/regolith_melt_zone.f90\
regolith/regolith_quadratic_func.f90\
regolith/regolith_quartic_func.f90\
regolith/regolith_mix.f90\
regolith/regolith_depth_model.f90\
regolith/regolith_streamtube.f90\
regolith/regolith_streamtube_head.f90\
regolith/regolith_streamtube_lineseg.f90\
regolith/regolith_subpixel_streamtube.f90\
regolith/regolith_transport.f90\
regolith/regolith_traverse_streamtube.f90\
regolith/regolith_subcrater_mix.f90\
regolith/regolith_melt_glass.f90\
regolith/regolith_superdomain.f90\
regolith/regolith_melt_zone_superdomain.f90\
regolith/regolith_streamtube_volume_func.f90\
regolith/regolith_shock_damage_zone.f90\
regolith/regolith_shock_damage.f90\
porosity/porosity_form_interior.f90\
realistic/realistic_perlin_noise.f90\
main/CTEM.f90
CLEANFILES = *.mod
