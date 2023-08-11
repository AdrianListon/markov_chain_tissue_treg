# write_dwell_time.r 

write_dwell_time <- function( celltype = celltype, ps, n.iter = 0,
                              sel_models_file_name, q12, states )

{
  pmc <- get_parabiosis_const()
  pft <- parabiosis_flow_titles( celltype = celltype )
  
  dAllTissues_Models_HDI <- get.dAllTissues_Models_HDI_simple( 
    celltype = celltype, selected.models.file = 
      sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name ) ) 
  
  model.id <- gsub( pattern = ".*_([0-9]{3,4}[a-z]{0,3})_.*", replacement = "\\1",
                    x = sel_models_file_name )
  
  
  badformat_dwell_csv_full <- sprintf( 
    "%s/%s/%s/%s.csv", ps$RESULTS_PATH, celltype, "Figures", sprintf( 
      "f%s_%s_%s_Q12=%s_%s_Dwell_Time", model.id, celltype, n.iter, q12, states ) )
  
  dAllTissues_Models_HDI %>% 
    left_join( pmc$dTissueAllOrderedGroup.f %>%                                  
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    # multi-tissue models
    mutate( par_orig = par ) %>%                                                 
    mutate( par = gsub( pattern = "(q[0-9]{2,2})(_i\\[1\\])", replacement = "\\1", 
                        x = par_orig ) ) %>% 
    filter( par %in% c( "q44", "q55", "q66" ) ) %>% 
    mutate( dwell_time_mode_DAYS = -1 / mode,
            dwell_time_hdi_80_upper = -1 / hdi.80.lower,
            dwell_time_hdi_80_lower = -1 / hdi.80.upper,
            dwell_time_hdi_50_upper = -1 / hdi.50.lower,
            dwell_time_hdi_50_lower = -1 / hdi.50.upper ) %>% 
    mutate( cellstate = ifelse( par == "q66", "CD69+",
                                ifelse( par == "q55", "Antigen-experienced",
                                        ifelse( par == "q44", "Naive", "OTHER" ) ) ) ) %>% 
    dplyr::select( celltype, tissue, cellstate, starts_with( "dwell_time" ), par ) %>% 
    arrange( desc( par ), tissue ) %>% 
    dplyr::distinct() %>% 
    write_csv( ., badformat_dwell_csv_full ) 
  
    format_dwell_time_table( csv_name = badformat_dwell_csv_full, ps = ps )
    plot_dwell_times( ps = ps, csv_name = badformat_dwell_csv_full, 
                           celltype = celltype, n.iter = n.iter, 
                           sel_models_file_name = sel_models_file_name,
                           q12 = q12, states = states )
}
