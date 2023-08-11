# plot_figures_2_4.r

plot_figures_2_4 <- function( celltype = celltype, ps, n.iter = 0, 
                              sel_models_file_name, q12, states )
{
  pmc <- get_parabiosis_const()
  pft <- parabiosis_flow_titles( celltype = celltype )
  
  dAllTissues_Models_HDI <- get.dAllTissues_Models_HDI_simple( 
    celltype = celltype, selected.models.file = 
      sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name ) )               
  
  model.id <- gsub( pattern = ".*_([0-9]{3,4}[a-z]{0,3})_.*", replacement = "\\1",
                    x = sel_models_file_name )
  
  dAllTissues_Models_HDI %>% 
    # For the vertical lines separating tissue groups:
    left_join( pmc$dTissueAllOrderedGroup.f %>%                                   
                 dplyr::select( tissue.all, f.tissue, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>% 
    # multi-tissue models
    mutate( par_orig = par ) %>%                                                 
    mutate( par = gsub( pattern = "(q[0-9]{2,2})(_i\\[1\\])", replacement = "\\1", 
                        x = par_orig ) ) %>% 
    
    left_join( ., pft$dFlowsGroups %>%
                 dplyr::select( par, f.flow.subgroup, flow_between ),
               by = c( "par" = "par" ) ) %>%
    # get titles   
    left_join( ., pft$dFlowsTitles %>% dplyr::select( flow_between, flow.title ),
               by = c( "flow_between" = "flow_between" ) ) %>%
    
    mutate( f.flow.subgroup.to.plot = ifelse( 
      flow_between %in% c( "1_to_2", "4_to_1", "5_to_2", "6_to_5" ), 
      "not.to.plot", f.flow.subgroup ) ) %>% 
    group_by( f.flow.subgroup.to.plot ) %>%                                      
    mutate( hdi.80.lower.flow.subgroup.min = min( hdi.80.lower, na.rm = TRUE ),
            hdi.80.upper.flow.subgroup.max = max( hdi.80.upper, na.rm = TRUE )
    ) %>%
    ungroup() -> 
    dTMP.2
  
  
  xintercept1 <- dTMP.2$f.tissue[ dTMP.2$f.tissue.group == "Lymphoid" ] %>% 
    unique() %>% length()
  xintercept2 <- dTMP.2$f.tissue[ dTMP.2$f.tissue.group == "Non-lymphoid" ] %>% 
    unique() %>% length()
  
  l_ggfigs <- list()
  FLOWS <- dTMP.2$flow_between %>% na.omit %>% unique %>% sort

  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  for ( i_flow in FLOWS ) {
    MIN_0 <- 0.0001
    hdi.80.lower.flow.subgroup.min.this <- 
      min( dTMP.2$hdi.80.lower.flow.subgroup.min[ 
        dTMP.2$flow_between == i_flow ], na.rm = TRUE )
    if( hdi.80.lower.flow.subgroup.min.this == 0 ) 
      hdi.80.lower.flow.subgroup.min.this <- MIN_0
    
    hdi.80.upper.flow.subgroup.max.this <- 
      max( dTMP.2$hdi.80.upper.flow.subgroup.max[ 
        dTMP.2$flow_between == i_flow ], na.rm = TRUE )
    if( hdi.80.upper.flow.subgroup.max.this == 0 ) 
      hdi.80.upper.flow.subgroup.max.this <- MIN_0
    
    dTMP.2 %>% 
      filter( !is.na( f.tissue ) ) %>% 
      filter( flow_between == i_flow ) %>%
      mutate( hdi.80.lower = ifelse ( hdi.80.lower < MIN_0, MIN_0, hdi.80.lower ),
              hdi.50.lower = ifelse ( hdi.50.lower < MIN_0, MIN_0, hdi.50.lower ) ) %>% 
      
      ggplot( aes( x = f.tissue, y = mode * 1e3 ) ) +                              
      geom_point( size = 3, alpha = 0.5 ) +
      scale_y_log10( n.breaks = 6,
                     limits = c( hdi.80.lower.flow.subgroup.min.this * 1e3 ,  
                                 hdi.80.upper.flow.subgroup.max.this * 1e3 ),    
                     labels = scales::trans_format( 
                       "log10", scales::math_format( 10^.x ) ) ) +               
      annotation_logticks( sides = "l" ) +
      geom_crossbar( aes( ymin = hdi.80.lower * 1e3, ymax = hdi.80.upper * 1e3 ), 
                     colour = "grey", width = 0.6 ) +
      geom_crossbar( aes( ymin = hdi.50.lower * 1e3, ymax = hdi.50.upper * 1e3 ), 
                     colour = "blue" ) +
      
      # graphics
      geom_vline( xintercept = c( xintercept1 + 0.5, xintercept1 + xintercept2 + 0.5 ) ) +
      geom_hline( yintercept = 1, colour = "darkred", linetype = 2, size = 1 )  +
      labs( title = dTMP.2$flow.title[ dTMP.2$flow_between == i_flow ] %>% 
              na.omit %>% unique(),
            y = "Rate (events/1000 cells/day)" ) +     
      theme_classic() +
      theme( plot.title = element_text( size = 12 * f, hjust = 0 ), 
             axis.title.x = element_blank(),
             axis.title.y = element_text( size = 12 * f, hjust = 1 ),  
             
             axis.text.x = element_text( angle = 45, hjust = 1, size = 12 * f ),
             axis.text.y = element_text( size = 12 * f ),  
             
             legend.text = element_text( size = 12 * f ),
             legend.title = element_text( size = 12 * f ),
             
             panel.border = element_blank(),             
             panel.grid.major = element_blank(),         
             panel.grid.minor = element_blank() ) ->
      ggfig; ggfig

    l_ggfigs[[ which( i_flow == FLOWS ) ]] <- ggfig
    
    fig.filename.start <- sprintf( 
      "f%s_%s_%s_Q12=%s_%s_Mode_HDI", model.id, celltype, n.iter, q12, states )
    
    path.figures <- sprintf( "%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures" )
    if ( !dir.exists( path.figures ) ) { 
      dir.create( path.figures, recursive = TRUE ) }
    path.fig.filename.start <- sprintf( "%s/%s", path.figures, fig.filename.start )
    if ( !dir.exists( path.fig.filename.start ) ) {
      dir.create( path.fig.filename.start, recursive = FALSE ) }
  }
  
  if ( FIG_5BCD <- TRUE ) {
    gg.fig.5BCD <- ggpubr::ggarrange( 
      l_ggfigs[[ which( "1_to_4" == FLOWS ) ]], 
      l_ggfigs[[ which( "2_to_5" == FLOWS ) ]], 
      l_ggfigs[[ which( "3_to_6" == FLOWS ) ]],
      labels = c( "B", "C", "D" ),
      ncol = 3, nrow = 1 )
    

    ggpubr::ggexport( gg.fig.5BCD,
                      filename = sprintf( "%s/%s_Figure_5BCD.pdf",
                                          path.figures, fig.filename.start ),
                      width = 21 * cm_to_in, height = 6 * cm_to_in )
    
    add_plot_to_figs_rda( gg_add = gg.fig.5BCD, gg_add_name = "gg.fig.5BCD",
                          ps = ps, celltype = celltype )
  }
  
  if ( FIG_S5B <- TRUE ) {
    gg.fig.S5B <- ggpubr::ggarrange( 
      l_ggfigs[[ which( "1_to_4" == FLOWS ) ]], 
      l_ggfigs[[ which( "2_to_5" == FLOWS ) ]], 
      l_ggfigs[[ which( "3_to_6" == FLOWS ) ]],
      
      l_ggfigs[[ which( "4_to_1" == FLOWS ) ]], 
      l_ggfigs[[ which( "5_to_2" == FLOWS ) ]], 
      l_ggfigs[[ which( "6_to_3" == FLOWS ) ]], 
      
      l_ggfigs[[ which( "4_to_1" == FLOWS ) ]], 
      l_ggfigs[[ which( "5_to_1" == FLOWS ) ]], 
      l_ggfigs[[ which( "6_to_1" == FLOWS ) ]],
      
      labels = c( LETTERS[ 1:9 ] ),
      ncol = 3, nrow = 3 )
    
    ggpubr::ggexport( gg.fig.S5B,
                      filename = sprintf( "%s/%s_Figure_S5B.pdf", 
                                          path.figures, fig.filename.start ), 
                      width = 12, height = 8 )
  }
}  
  
