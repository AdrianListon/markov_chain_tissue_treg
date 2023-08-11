# load_my_functions.r

load_my_functions <- function( ps )
{
  source( sprintf( "%s/%s", ps$CODE_PATH, "get_parabiosis_const.r" ) )
  
  source( sprintf( "%s/%s", ps$CODE_PATH, "estimate_one_markov_model.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "prep_for_markov_model.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "get_HDI.r" ) )               
  source( sprintf( "%s/%s", ps$CODE_PATH, "calculate_HDI.r" ) )               
  source( sprintf( "%s/%s", ps$CODE_PATH, "markov_premodel_prints_eps.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "select_mcmc_chains.r" ) )
  
  source( sprintf( "%s/%s", ps$CODE_PATH, "add_plot_to_figs_rda.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "figure_flow_titles.r" ) )               
  source( sprintf( "%s/%s", ps$CODE_PATH, "preproc_fs_for_figures_2_4_3_5.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "plot_figures_2_4.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "group_flows.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "write_dwell_time.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "format_dwell_time_table.r" ) )
  source( sprintf( "%s/%s", PS$CODE_PATH, "visualize_results.r" ) )  
  
  source( sprintf( "%s/%s", ps$CODE_PATH, "converting_Q_and_flow_tables.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "prep_data_figures_6.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "merge_fig5_subfigures.r" ) )
  source( sprintf( "%s/%s", ps$CODE_PATH, "plot_all_tissue_trajectories.r" ) )
}