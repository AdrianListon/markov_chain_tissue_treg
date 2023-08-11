# visualize_results.r

visualize_results <- function( 
    celltype, ps, n.iter, sel_models_file_name, q12, states,
    model.name, model.ver, setup.mcmc.fname )
{
  plot_figures_2_4( 
    celltype = celltype, ps = ps, n.iter = n.iter,
    sel_models_file_name = sel_models_file_name,
    q12 = q12, states = states )
  # hack: ceiling(...) instead of round() which is not OS compatible    
  
  write_dwell_time( 
    celltype = celltype, ps = ps, n.iter = n.iter,
    sel_models_file_name = sel_models_file_name,
    q12 = q12, states = states )
  
  plot_aggr_tissuegroup_cellstate_barplot( 
    celltype = celltype, ps = ps )
  
  write.aggregated.group.flows( 
    celltype = celltype, ps = ps, n.iter = n.iter,
    sel_models_file_name = sel_models_file_name,
    q12 = q12, states = states,
    model.name = model.name, model.ver = model.ver, 
    setup.mcmc.fname = setup.mcmc.fname, write_csv = TRUE )
  
  get_aggr_tissuegroup_cellstate_areas( 
    celltype = celltype, write_csv = TRUE )
  
  get_aggr_total_counts_areas( 
    celltype = celltype, write_csv = TRUE )
  
  merge_fig5_subfigures( 
    celltype = celltype, ps = ps, n.iter = n.iter,
    sel_models_file_name = sel_models_file_name,
    q12 = q12, states = states )
  
  plot_all_tissue_trajectories(
    celltype = celltype, ps = ps, n.iter = n.iter,
    sel_models_file_name = sel_models_file_name,
    q12 = q12, states = states )
}