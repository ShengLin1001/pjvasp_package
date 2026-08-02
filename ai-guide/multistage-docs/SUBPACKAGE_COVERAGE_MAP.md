# 子包覆盖地图


## build（11 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `build/bulk/create.py` | 359 | 11 | 0 | create_fcc_111, create_hcp_basal, create_hcp_prism1, create_hcp_prism2, create_supercell_fast, vasp_create_hcp_basal |
| `build/bulk/gsfe.py` | 77 | 1 | 0 | create_gsfe_model |
| `build/film/extrfilm.py` | 90 | 4 | 0 | my_extr_lattice, my_extr_thick, my_extr_etot, cal_area |
| `build/film/findcubic.py` | 201 | 3 | 0 | find_cubic, find_optimal_cell_shape, get_deviation_from_optimal_cell_shape |
| `build/film/findhetero.py` | 329 | 8 | 0 | find_hetero, my_plot_results, set_length, magnitude, build_supercells, stack_atoms |
| `build/film/findprim.py` | 169 | 2 | 0 | my_find_prim, check_direction |
| `build/film/hydroxyl.py` | 457 | 5 | 0 | add_hydro_atoms, get_neighbors, find_matching_atom_in_bulk, find_unsaturated_neighbors, passivate_surface_custom |
| `build/film/stretch.py` | 515 | 8 | 0 | generate_film, stretch_list_along_direction_to_cell, stretch_along_direction_to_cell, print_after_what, file_name, my_find_num_per_slab |
| `build/workflow/general.py` | 54 | 3 | 0 | cp_vaspfiles, rm_i, compare_three_lattices |
| `build/workflow/hoec.py` | 390 | 3 | 0 | prepare_hoec_reference, deform_atoms, generate_hoec_dirs |
| `build/workflow/kpar_ncore.py` | 336 | 11 | 0 | get_default_pairs, check_pairs, check_source_inputs, copy_case_inputs, write_case_control_files, change_case_incar |

## calculate（12 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `calculate/calenergy/surfenergy.py` | 77 | 1 | 0 | cal_surface_energy |
| `calculate/calmath/matrix.py` | 68 | 1 | 0 | hermite_normal_form |
| `calculate/calmechanics/deformation.py` | 32 | 1 | 0 | cal_deform_matrix |
| `calculate/calmechanics/hoec.py` | 1186 | 24 | 1 | rot_z, rot_axis, voigt_strain_M, close_group, monomials, mono_eval_vec |
| `calculate/calmechanics/strain.py` | 219 | 5 | 0 | cal_principal_and_shear_strain_root, cal_strain_matrix_root, cal_strain_matrix, cal_principal_and_shear_strain, cal_von_mises_strain |
| `calculate/calmechanics/stress.py` | 1 | 0 | 0 |  |
| `calculate/calmechanics/stretch.py` | 92 | 2 | 0 | cal_relative_stretch, cal_stretch |
| `calculate/calmismatch/calhetero.py` | 271 | 6 | 0 | compare_atoms, cal_mismatch, cal_stretch_lattice, relative_diff, cal_atom_num, filter_results |
| `calculate/calqm/kpoints.py` | 206 | 6 | 0 | get_kpoints_by_size, get_size_by_distance, get_distance_by_size, cal_reciprocal_matrix, cal_reciprocal_matrix2, get_lattice_information |
| `calculate/electronic_structure/plotter.py` | 4586 | 10 | 6 | plot_fermi_surface, plot_wigner_seitz, plot_lattice_vectors, plot_path, plot_labels, fold_point |
| `calculate/electronic_structure/universal.py` | 155 | 3 | 0 | get_n_band, summarize_band_structure_info, get_n_kpoints_band |
| `calculate/material_science/schmid.py` | 80 | 1 | 0 | cal_fcc_schmid_factors |

## io（5 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `io/extxyz.py` | 49 | 1 | 0 | extxyz_to_atomlist |
| `io/general.py` | 149 | 3 | 0 | general_read, general_write, get_col_type |
| `io/post/construct.py` | 36 | 1 | 0 | create_content |
| `io/post/write.py` | 29 | 1 | 0 | write_content_to_file |
| `io/vasp.py` | 445 | 8 | 0 | my_read_vasp, my_write_vasp, read_vasp_list, get_atomtypes, atomtypes_outpot, get_atomtypes_from_formula |

## post（12 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `post/Cij_energy.py` | 324 | 7 | 0 | post_lammps_Cij_energy, check_ldata, fit_cij_energy, write_cij_energy, read_cij_energy, read_deform_data |
| `post/convergence.py` | 170 | 3 | 0 | post_convergence, my_write_convergence, my_read_convergence |
| `post/E_in_1_2_bulk.py` | 133 | 3 | 0 | post_E_in_1_2_bulk, my_write_E_in_1_2_bulk, my_read_E_in_1_2_bulk |
| `post/general.py` | 148 | 4 | 0 | my_sort, my_ployfit, my_read_y_dir_contcar, get_structure_info |
| `post/gsfe.py` | 432 | 5 | 0 | post_gsfe, check_constraints, find_sf_usf, write_output, read_output |
| `post/hoec_energy.py` | 1050 | 17 | 0 | run_univ_post, load_hoec_manifest, read_element, collect_mode_data, fit_P, fit_P_fixed |
| `post/kpar_ncore.py` | 365 | 9 | 0 | run_univ_post, read_kpar_ncore_times, read_kpar_ncore_energies, read_vasp_natoms, get_kpar_ncore_natoms, get_delta_energies |
| `post/neb.py` | 490 | 12 | 0 | post_neb, my_copy_neb_files, my_write_neb, my_add_head, my_read_neb, my_spline_neb |
| `post/newmain.py` | 795 | 7 | 9 | post_general, myfindall, read_OUTCAR, check_convergence, write_convergence, check_same |
| `post/oldmain.py` | 1285 | 29 | 9 | myfindall, my_extr_lattice, my_extr_thick, my_extr_etot, read_OUTCAR, read_vasp_list |
| `post/relax_convergence.py` | 139 | 4 | 0 | my_univ_post_convergence, read_univ_post_convergence, _to_float, _to_int |
| `post/stretch.py` | 304 | 5 | 0 | post_stretch, post_lammps_stretch, get_stretch_type, my_write_stretch, my_read_stretch |

## universal（28 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `universal/atom/delatom.py` | 63 | 2 | 0 | mydel_pos_type, check_position |
| `universal/atom/density.py` | 18 | 1 | 0 | cal_density |
| `universal/atom/fixatom.py` | 38 | 1 | 0 | fixatoms |
| `universal/atom/miller.py` | 78 | 1 | 0 | three_index_to_four_index |
| `universal/atom/moveatom.py` | 108 | 1 | 0 | move_atoms |
| `universal/atom/neighbor.py` | 32 | 1 | 0 | get_neighbor_distances |
| `universal/atom/rotate.py` | 32 | 1 | 0 | cal_rotate |
| `universal/check/atoms.py` | 175 | 4 | 0 | get_cna_count, get_cna_count_dict, get_cna_count_vasp, check_phase_transition |
| `universal/check/checkinput.py` | 36 | 1 | 0 | check_input |
| `universal/check/convergence.py` | 59 | 2 | 0 | check_and_submit_jobs_in_ydir, check_and_submit_jobs_in_subdir |
| `universal/check/type.py` | 28 | 4 | 0 | check_positive_int, check_none, check_basename, check_absolute_path |
| `universal/data/dataadjust.py` | 95 | 5 | 0 | rm_blank, my_add_list, normalize_float_int, myjust, my_down_up |
| `universal/data/datachange.py` | 93 | 3 | 0 | list_to_char, char_to_dic, dic_to_char |
| `universal/data/patterntrans.py` | 42 | 1 | 0 | my_pattern_trans |
| `universal/math/operations.py` | 60 | 2 | 0 | get_integration, get_integration_curve |
| `universal/matrix/adjust.py` | 39 | 1 | 0 | adjust_matrix |
| `universal/plot/atominfo.py` | 128 | 3 | 0 | my_plot_interlayer_distance, my_plot_zpositions, my_plot_rdf |
| `universal/plot/energy.py` | 230 | 1 | 0 | my_plot_energy_components |
| `universal/plot/general.py` | 1203 | 19 | 0 | check_font_size, check_axes_size, check_all_rcParams, get_ploted_figure, get_points_on_markers_boundary, add_color_band |
| `universal/plot/n2p2.py` | 668 | 19 | 0 | _tag_colors, _format_plain, my_plot_learning_curve, my_plot_compare, _bar_yerr, my_plot_rmse_by_tag |
| `universal/plot/oldplotdos.py` | 302 | 5 | 0 | my_plot_horizontal_vertical, my_plot_orientation, my_plot_complete_dos, my_plot_idos, my_plot_element_spd_dos |
| `universal/plot/plot.py` | 471 | 5 | 0 | general_modify_ploted_figure, my_plot, my_plot_brokenaxed, my_plot_modify_ploted_figure, my_plot_colorbar |
| `universal/plot/plotting.py` | 742 | 11 | 0 | pretty_plot, pretty_plot_two_axis, pretty_polyfit_plot, _decide_fontcolor, periodic_table_heatmap, format_formula |
| `universal/plot/render.py` | 56 | 1 | 0 | my_render |
| `universal/plot/workflow.py` | 1222 | 14 | 0 | my_plot_convergence, my_plot_cij_energy, my_plot_gsfe, my_plot_gsfe_displacement, _apply_relax_yscale, my_plot_relax_convergence |
| `universal/print/print.py` | 52 | 5 | 0 | pr, er, warn, fail, confirm_prepare_outdir |
| `universal/print/printafter.py` | 63 | 4 | 0 | print_after_read, print_after_blank, print_after_not_supported, print_after_cant_read |
| `universal/search/find.py` | 56 | 2 | 0 | find_line_position, extract_line_at_position |

## ml（13 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `ml/confusionmatrix.py` | 97 | 0 | 1 | ConfusionMatrix |
| `ml/dataset.py` | 180 | 2 | 1 | get_all_data, get_mean_std, CustomDataset |
| `ml/model.py` | 140 | 3 | 0 | train_model, visualize_model, visualize_model_predictions |
| `ml/n2p2/calculate/cur.py` | 211 | 4 | 0 | _read_function_data, collect_sf_features, filter_zero_columns, cur_select |
| `ml/n2p2/calculate/post.py` | 152 | 6 | 0 | read_learning_curve, read_normalization, read_trainpoints, read_trainforces, rmse_me, build_rmse_by_tag_df |
| `ml/n2p2/calculate/sf.py` | 890 | 9 | 1 | load_from_input_nn, get_largest_rc_from_input_nn, get_radial_pairs, get_angular_pairs, get_leta_lshift_from_N, _format_sf_block |
| `ml/n2p2/dataset.py` | 660 | 4 | 1 | _choose_inplane_rep, _is_orthogonal_cell, _stretch_lattice_values_for_epoch_convention, read_dft_reference, nnpdata |
| `ml/n2p2/sfparamgen/generate-acsf-angular_narrow.py` | 44 | 0 | 0 |  |
| `ml/n2p2/sfparamgen/generate-acsf-angular_wide.py` | 45 | 0 | 0 |  |
| `ml/n2p2/sfparamgen/generate-acsf-radial.py` | 39 | 0 | 0 |  |
| `ml/n2p2/sfparamgen/sfparamgen.py` | 724 | 0 | 1 | SymFuncParamGenerator |
| `ml/n2p2/workflow.py` | 1972 | 0 | 1 | PeiN2p2 |
| `ml/plot.py` | 25 | 1 | 0 | my_imshow |

## slurm（1 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `slurm/submit.py` | 874 | 13 | 0 | _write_slurm_script, check_wall_time, generate_script_header, generate_launcher_command, generate_slurm_script_base, generate_slurm_script_sequential |

## cr（2 个模块）

| 文件 | LOC | 函数数 | 类数 | 公开函数/类 |
|------|-----|--------|------|-------------|
| `cr/crplotkcontactgraph.py` | 359 | 7 | 2 | get_path, get_data, H_handle, Freq_handle, K_handle, E_handle |
| `cr/crsum3.py` | 916 | 10 | 1 | func, get_average_K_Path, draw_point, handle_data, get_data_sum_K, plot |