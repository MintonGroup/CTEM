bin_PROGRAMS = CTEM 
#ifort optimized flags
#AM_FCFLAGS =  -O3 -openmp -parallel -xHost -ipo -assume byterecl -m64 -heap-arrays -FR
#ifort debug flags
#AM_FCFLAGS = -O0 -p -g -openmp -debug all -traceback -CB

#gfortran optimized flags
AM_FCFLAGS =  -O3 -fopenmp -ffree-form
#gfortran debug flags
#AM_FCFLAGS =  -O0 -g -fopenmp -fbounds-check -Wall -Warray-bounds -Warray-temporaries -Wimplicit-interface -ffree-form

CTEM_SOURCES = globals/module_globals.f90\
util/module_util.f90\
io/module_io.f90\
crater/module_crater.f90\
ejecta/module_ejecta.f90\
init/module_init.f90\
seismic/module_seismic.f90\
regolith/module_regolith.f90\
porosity/module_porosity.f90\
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
io/io_read_const.f90\
io/io_get_token.f90\
io/io_input.f90\
io/io_read_prod.f90\
io/io_read_vdist.f90\
io/io_read_surf.f90\
io/io_read_regotrack.f90\
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
ejecta/ejecta_emplace.f90\
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
crater/crater_find_visible.f90\
crater/crater_averages.f90\
crater/crater_emplace.f90\
crater/crater_form_interior.f90\
crater/crater_form_exterior.f90\
crater/crater_record.f90\
crater/crater_tally_true.f90\
crater/crater_tally_observed.f90\
crater/crater_tally_calibrated_count.f90\
crater/crater_slope_collapse.f90\
crater/crater_soften.f90\
crater/crater_subpixel_diffusion.f90\
init/init_domain.f90\
init/init_dist.f90\
init/init_surf.f90\
init/init_regolith_stack.f90\
init/init_regolith_parab.f90\
seismic/seismic_shake.f90\
seismic/seismic_distance.f90\
seismic/seismic_table_define.f90\
seismic/seismic_kdiff_func.f90\
regolith/regolith_circle_sector_func.f90\
regolith/regolith_cubic_func.f90\
regolith/regolith_melt_func.f90\
regolith/regolith_melt_fraction.f90\
regolith/regolith_melt_zone.f90\
regolith/regolith_pop.f90\
regolith/regolith_push.f90\
regolith/regolith_quadratic_func.f90\
regolith/regolith_quartic_func.f90\
regolith/regolith_mix.f90\
regolith/regolith_mix_porous_regime.f90\
regolith/regolith_streamtube.f90\
regolith/regolith_streamtube_head.f90\
regolith/regolith_streamtube_lineseg.f90\
regolith/regolith_subpixel_streamtube.f90\
regolith/regolith_transport.f90\
regolith/regolith_traverse_pop.f90\
regolith/regolith_traverse_streamtube.f90\
porosity/porosity_form_interior.f90\
main/CTEM.f90
CLEANFILES = *.mod
