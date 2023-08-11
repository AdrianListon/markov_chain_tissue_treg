# format_dwell_time_table.r

format_dwell_time_table <- function( csv_name, ps = NA )
{
  pmc <- get_parabiosis_const()  
  
  read_csv( csv_name ) %>% 
    dplyr::select( -c( dwell_time_hdi_80_lower, dwell_time_hdi_80_upper, 
                       dwell_time_hdi_50_lower, dwell_time_hdi_50_upper, 
                       par ) ) %>% 
    distinct() %>% 
    mutate( dwell_time_mode_DAYS = round( dwell_time_mode_DAYS, digits = 1 ) ) %>% 
    left_join( ., pmc$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    tidyr::spread( cellstate, dwell_time_mode_DAYS ) %>% 
    dplyr::rename( cd69p = `CD69+`, activ = `Antigen-experienced`, naive = Naive ) %>% 
    arrange( f.tissue.group ) %>% 
    mutate( tissue_group = f.tissue.group %>% as.character ) %>% 
    dplyr::select( celltype, tissue, tissue_group, cd69p, activ, naive, 
                   f.tissue.group ) ->
    d_tissues
    
  d_tissues %>% 
    group_by( f.tissue.group ) %>% 
    summarise( celltype = first( celltype ),
               tissue = "Mean",
               cd69p = mean( cd69p ) %>% round( digits = 1 ),
               activ = mean( activ ) %>% round( digits = 1 ), 
               naive = mean( naive ) %>% round( digits = 1 ) ) %>% 
    mutate( tissue_group = f.tissue.group %>% as.character ) %>% 
    dplyr::select( celltype, tissue, tissue_group, cd69p, activ, naive ) ->
    d_tissue_group_means
    
  d_tissues %>% 
    dplyr::select( celltype, tissue, tissue_group, cd69p, activ, naive ) %>% 
    rbind( ., "", "", d_tissue_group_means ) %>% 
    write_csv( file = sprintf( "%s_NICE.csv", sub( '\\.csv$', '', csv_name ) ) )
}


plot_dwell_times <- function( csv_name, ps = NA, celltype, 
                              plot = TRUE, plot_CrI = TRUE,
                              n.iter, sel_models_file_name, q12, states ) 
{
  pmc <- get_parabiosis_const() 
  
  model.id <- gsub( pattern = ".*_([0-9]{3,4}[a-z]{0,3})_.*", replacement = "\\1",
                    x = sel_models_file_name )
  
  read_csv( csv_name ) %>% 
    dplyr::select( -par ) %>% 
    distinct() %>% 
    left_join( ., pmc$dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    # rename for Tissue Treg paper:
    mutate( cellstate = ifelse( 
      cellstate == "Naive", "resting", cellstate ) ) %>% 
    mutate( cellstate = ifelse( 
      cellstate == "Antigen-experienced", "activated", cellstate ) ) %>% 
    mutate( cellstate = factor(
      cellstate, levels = c( "resting", "activated", "CD69+" ) ) ) ->
    dDwellTime
  
  xintercept1 <- dDwellTime$f.tissue[ dDwellTime$f.tissue.group == "Lymphoid" ] %>% 
    unique() %>% length()
  xintercept2 <- dDwellTime$f.tissue[ dDwellTime$f.tissue.group == "Non-lymphoid" ] %>% 
    unique() %>% length()
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  dDwellTime %>% 
    ggplot( aes( x = f.tissue, y = dwell_time_mode_DAYS, fill = cellstate ) ) +
    geom_bar( stat = "identity", position = "dodge" ) +
    geom_errorbar( aes( ymin = dwell_time_hdi_80_upper, 
                        ymax = dwell_time_hdi_80_lower ), width = 0.2,
                        position = position_dodge( 0.9 ) ) + 
    geom_vline( xintercept = c( xintercept1 + 0.5, 
                                xintercept1 + xintercept2 + 0.5 ),
                linetype = "dashed" ) +
    labs( y = "Mean Dwell Time (days)", fill = "Cell state" ) +    
    theme_classic() +
    theme( plot.title = element_text( size = 12 * f, hjust = 0 ), 
           axis.title.x = element_blank(),
           axis.title.y = element_text( size = 12 * f, hjust = 0.5 ),  
           axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),
           axis.text.y = element_text( size = 12 * f ),  
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.text = element_text( size = 12 * f ),
           legend.title = element_text( size = 12 * f ),
           strip.text.x = element_text( size = 12 * f ) ) -> 
  gg.fig.5E
  
  if ( plot ) {
    fig.filename.start <- sprintf( 
      "f%s_%s_%s_Q12=%s_%s_Mode_HDI", model.id, celltype, n.iter, q12, states )
    
    ggpubr::ggexport( gg.fig.5E, filename = sprintf( 
      "%s/%s/%s_Fig_5E_dwelltimes.pdf", ps$RESULTS_PATH, celltype, 
      fig.filename.start ), width = 21 * cm_to_in, height = 6 * cm_to_in ) 
    add_plot_to_figs_rda( gg_add = gg.fig.5E, gg_add_name = "gg.fig.5E",
                          ps = ps, celltype = celltype )
  }
}
