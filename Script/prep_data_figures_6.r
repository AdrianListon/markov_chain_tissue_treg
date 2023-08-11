# prep_data_figures_6.r                               


write.aggregated.group.flows <- function( 
    ps = PS, celltype, write_csv = FALSE,
    n.iter, sel_models_file_name, q12, states,
    model.name, model.ver, setup.mcmc.fname )

{
  model.id <- gsub( pattern = ".*_([0-9]{3,4}[a-z]{0,3})_.*", replacement = "\\1",
                    x = sel_models_file_name )
  
  get.tissues.average.flows( ps = ps, celltype = celltype, flow.part = "exit",
                             model.name = model.name, model.ver = model.ver, 
                             setup.mcmc.fname = setup.mcmc.fname ) %>%     
    bind_rows( ., get.tissues.average.flows( 
      ps = ps, celltype = celltype, flow.part = "entry",
      model.name = model.name, model.ver = model.ver, 
      setup.mcmc.fname = setup.mcmc.fname ) ) ->
    average.flows.celltype
  
  
  average.flows.celltype %>%
    mutate( arrow.group = sprintf( 
      "g_%s_%s", 
      ifelse( population.1 < population.2, population.1, population.2 ),
      ifelse( population.1 > population.2, population.1, population.2 ) ) ) %>% 
    
    mutate( median.flow.1000.round = round( 1000 * median.flow, 0 ) ) %>% 
    dplyr::select( f.tissue.group, arrow.group, population.1, population.2, par, 
                   mean.flow, flow.mean.w.by.tissue.in.tissue.group, flow.part, 
                   median.flow, median.flow.1000.round, log10.median.flow, log10.median.flow.1000 ) %>% 
    
    mutate( log10.median.flow.1000.tmp.3 = 
              ifelse( log10.median.flow.1000 < 1, 1, 
                      ifelse( log10.median.flow.1000 > 4, 
                              4, log10.median.flow.1000 ) ) - 1 ) %>% 
    mutate( log10.median.flow.1000.deg.3 = 
              log10.median.flow.1000.tmp.3 / 3 * 90 ) %>% 
    
    mutate( log10.median.flow.1000.tmp.4 = 
              ifelse( log10.median.flow.1000 < 0, 0, 
                      ifelse( log10.median.flow.1000 > 4, 
                              4, log10.median.flow.1000 ) ) ) %>% 
    mutate( log10.median.flow.1000.deg.4 = 
              log10.median.flow.1000.tmp.4 / 4 * 90 ) %>% 
    
    dplyr::select( -c( 
      log10.median.flow.1000.tmp.3, log10.median.flow.1000.tmp.4,
      log10.median.flow.1000.deg.4 ) ) %>% 
    
    
    arrange( arrow.group, f.tissue.group, par ) -> 
    aggregated.group.flows
  
  if ( write_csv ) {
    fig.filename.start <- sprintf( 
      "f%s_%s_%s_Q12=%s_%s_Mode_HDI", model.id, celltype, n.iter, q12, states )
    
    path.flow.figure <- sprintf( "%s/%s/Flow_Diagram", ps$RESULTS_PATH, celltype )
    if( !dir.exists( path.flow.figure ) ) { dir.create( path.flow.figure ) }
    write_csv( aggregated.group.flows, sprintf( 
      "%s/%s_flow_diagram_medians_means_degs.csv", 
      path.flow.figure, fig.filename.start ) )
    }
  return( aggregated.group.flows )
}


get_cellstate_areas <- function( 
    ps = PS, f_tissue = "Brain", celltype = "Treg", 
    f_week = 0, hd  = "host" ) 
{
  stopifnot( hd %in% c( "host", "donor" ) )
  
  parabio.file.name <- 
    sprintf( "%s/%s/%s_parabiosis_data_%s.csv", 
    ps$PROCESSED_PATH, celltype, celltype, 
    tolower( ifelse( f_tissue == "Blood", "Brain", f_tissue ) ) )
  parabio.data <- read.csv( parabio.file.name )
  colnames( parabio.data ) <- 
    gsub( pattern = "\\.(.*)\\.", 
          replacement = ".celltype.", colnames( parabio.data ) )
  if ( hd == "host" ) {
    parabio.data$sel.hd.celltype.naive <- parabio.data$host.celltype.naive
    parabio.data$sel.hd.celltype.activ <- parabio.data$host.celltype.activ
    parabio.data$sel.hd.celltype.cd69p <- parabio.data$host.celltype.cd69p
  } else {
    parabio.data$sel.hd.celltype.naive <- parabio.data$donor.celltype.naive
    parabio.data$sel.hd.celltype.activ <- parabio.data$donor.celltype.activ
    parabio.data$sel.hd.celltype.cd69p <- parabio.data$donor.celltype.cd69p
  }
    
  parabio.data %>% 
    filter( tissue == f_tissue, week == f_week ) %>% 
    summarise( area_naive = mean( sel.hd.celltype.naive ), 
               area_activ = mean( sel.hd.celltype.activ ), 
               area_cd69p = mean( sel.hd.celltype.cd69p ) ) %>%  
    sweep( ., 1, rowSums( . ), "/" ) ->
    areas_cellstates
  
  tibble( tissue = f_tissue, week = f_week, hd = hd,
          cellstate = sub( "area_", "", names( areas_cellstates ) ),
          area = as.numeric( areas_cellstates ),
          radius = sqrt( area ) ) ->
    dCellstateAreas

  return( dCellstateAreas )
}


get_aggr_tissuegroup_cellstate_areas <- function( 
    ps = PS, celltype, write_csv = FALSE )
{
  dTissueAllOrderedGroup.f <- get_parabiosis_const()$dTissueAllOrderedGroup.f
  
  d.tissue.all <- NULL
  for( tissue in dTissueAllOrderedGroup.f$tissue.all ) {
    d.tissue.all %>% 
      bind_rows( ., get_cellstate_areas( ps, tissue, celltype ) ) ->
      d.tissue.all
  }
  
  d.tissue.all %>% 
    left_join( ., dTissueAllOrderedGroup.f %>% 
                 dplyr::select( tissue.all, tissue.group ), 
               by = c( "tissue" = "tissue.all") ) %>% 
    group_by( tissue.group, cellstate ) %>% 
    summarise( median_tissue_group_area = median( area ), .groups = "drop" ) %>% 
    mutate( radius = sqrt( median_tissue_group_area ) ) %>% 
    # add Blood areas
    bind_rows( ., get_cellstate_areas( ps, "Blood", celltype ) %>%              
        dplyr::rename( median_tissue_group_area = area,
                       tissue.group = tissue ) ) %>% 
    arrange( tissue.group, cellstate ) ->
    median_tissuegroup_cellstate_areas
  
  if( write_csv ) { 
    path.flow.figure <- sprintf( "%s/%s/Flow_Diagram", ps$RESULTS_PATH, celltype )
    if( !dir.exists( path.flow.figure ) ) { dir.create( path.flow.figure ) }
    write_csv( median_tissuegroup_cellstate_areas, sprintf( 
      "%s/%s_tissuegroup_cellstates_areas.csv", path.flow.figure, celltype ) )
  }
  
  return( median_tissuegroup_cellstate_areas )
}


get_aggr_total_counts_areas <- function( ps = PS, celltype, write_csv = FALSE )
{
  total_counts <- read_csv( 
    sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv", 
             ps$PROCESSED_PATH, celltype ) ) 
  total_counts %>% 
    filter( Tissue %in% get_parabiosis_const()$tissue.all.ordered ) %>% 
    left_join( ., get_parabiosis_const()$dTissueAllOrderedGroupBlood.f %>% 
                 dplyr::select( tissue.all, tissue.group ), 
               by = c( "Tissue" = "tissue.all" ) ) %>% 
    group_by( tissue.group ) %>% 
    summarise( median_total_count = median( Mean ),
               mean_total_count = mean( Mean ) ) %>% 
    mutate( max_median_total_count = max( median_total_count ),
            rel_to_max_total_count = 
              median_total_count / max_median_total_count,
            radius_rel_total_count = sqrt( rel_to_max_total_count ) ) ->
    dTotalCountsAreas
  
  if( write_csv ) {   
    path.flow.figure <- sprintf( "%s/%s/Flow_Diagram", ps$RESULTS_PATH, celltype )
    if( !dir.exists( path.flow.figure ) ) { dir.create( path.flow.figure ) }
    write_csv( dTotalCountsAreas, sprintf( 
      "%s/%s_total_count_areas.csv", path.flow.figure, celltype ) )
  }

  return( dTotalCountsAreas )
}


plot_aggr_tissuegroup_cellstate_barplot <- function( 
  ps = PS, celltype = "Treg", plot = TRUE )
{

dCellstatesProps <- NULL
for ( tissue in get_parabiosis_const()$tissue.all.ordered )
  for ( hd in c( "host", "donor" ) )
    for ( week in c( 0, 1, 2, 4, 8, 12 ) ) {
      dCellstatesProps %>% 
        bind_rows( ., get_cellstate_areas( 
          ps = PS, f_tissue = tissue, celltype = "Treg", 
          hd = hd, f_week = week ) ) -> 
            dCellstatesProps
    }

dCellstatesProps %>% 
  filter( !( week == 0 & hd == "donor" ) ) %>%
  left_join( ., get_parabiosis_const()$dTissueAllOrderedGroupBlood.f,
             by = c( "tissue" = "tissue.all" ) ) %>% 
  group_by( tissue.group, hd, week, cellstate ) %>% 
  summarise( mean_prop = mean( area ) ) %>% 
  ungroup() %>% 
  bind_rows( ., tibble( 
    mean_prop = 1, hd = "donor", week = 0, 
    tissue.group = c( "Blood", "Non-lymphoid", "Lymphoid", "GALT" ) ) ) %>% 
  mutate( hd_wk = sprintf( "%s_wk%s", hd, str_pad( week, 2, pad = "0" ) ) ) %>%
  mutate( hd_ordered = ifelse( hd == "host", "01_host", "02_donor" ) ) %>%
  arrange( hd_ordered, week ) %>% 
  
  mutate( hd_hack = ifelse( hd == "host", "", " " ) ) %>% 
  mutate( hd_wk_ordered = factor( sprintf(
    "%s_wk%s", hd_ordered, str_pad( week, 2, pad = "0" ) ),
    labels = unique( sprintf( "%s%s%s", hd_hack, week, hd_hack ) ) ) ) %>%

  mutate( cellstate = ifelse( 
    cellstate == "naive", "resting", cellstate ) ) %>% 
  mutate( cellstate = ifelse( 
    cellstate == "activ", "activated", cellstate ) ) %>% 
  mutate( cellstate = ifelse( 
    cellstate == "cd69p", "CD69+", cellstate ) ) %>% 
  mutate( cellstate = factor(
    cellstate, levels = c( "resting", "activated", "CD69+" ) ) ) %>% 
  
  mutate( tissue.group = factor( 
    tissue.group, levels = c( "Blood", "Lymphoid", "Non-lymphoid", "GALT" ) ) ) ->
  dCellstatesPropsGroups
  
  cm_to_in <- 1 / 2.54   
  f <- cm_to_in * 1.7
  
  dCellstatesPropsGroups %>% 
  ggplot( ., aes( x = hd_wk_ordered, y = mean_prop, fill = cellstate ) ) + 
  geom_bar( position = "fill", stat = "identity" ) +
  facet_grid( cols = vars( tissue.group ) ) +
  theme_classic() +
  theme( plot.title = element_text( hjust = 0 ), 
         axis.title.x = element_blank(),
         axis.title.y = element_text( size = 12 * f, hjust = 0.5 ),  
         
         axis.text.x = element_text( size = 12 * f ),
         axis.text.y = element_text( size = 12 * f ),  
         
         panel.border = element_blank(),             
         panel.grid.major = element_blank(),         
         panel.grid.minor = element_blank(),
    
         legend.text = element_text( size = 12 * f ),
         legend.title = element_text( size = 12 * f ),
         strip.text.x = element_text( size = 12 * f, vjust = 0 ),
         strip.background = element_blank(),
         plot.tag = element_text( size = 12 * f, hjust = 0 ),
         plot.tag.position = c( 0.00, 0.037 ),
         plot.margin = margin( 5.5, 5.5, 11.5, 5.5, unit = "pt" ) ) +   
    
  scale_fill_discrete( na.value = "grey" ) + 
  
  annotate( "text", x = 3.5, y = -0.20, label = "host", size = 4 * f ) +
  annotate( "text", x = 9.5, y = -0.20, label = "donor", size = 4 * f ) +
  coord_cartesian( ylim = c( 0, NA ), clip = "off" ) +
  labs( y = "Cell state proportion", fill = "Cell state", tag = "week" ) -> 
  ggfig

  
  if ( plot ) {
    ggpubr::ggexport( ggfig, filename = sprintf( 
      "%s/%s/Fig_5A_cellstates_barchart.pdf", ps$RESULTS_PATH, celltype ), 
      width = 21 * cm_to_in, height = 6 * cm_to_in )
    
    add_plot_to_figs_rda( gg_add = ggfig, gg_add_name = "gg.fig.5A",
                          ps = ps, celltype = celltype )
    
  }
  return( ggfig )
}

