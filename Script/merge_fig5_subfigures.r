# merge_fig5_subfigures.r 

merge_fig5_subfigures <- function( 
    celltype, ps, n.iter = 0, sel_models_file_name, q12, states )

{
  pmc <- get_parabiosis_const()
  pft <- parabiosis_flow_titles( celltype = celltype )
  
  model.id <- gsub( pattern = ".*_([0-9]{3,4}[a-z]{0,3})_.*", replacement = "\\1",
                    x = sel_models_file_name )
  path.figures <- sprintf( "%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures" )
  if ( !dir.exists( path.figures ) ) { 
    dir.create( path.figures, recursive = TRUE ) }
  fig.filename.start <- sprintf( 
    "f%s_%s_%s_Q12=%s_%s", model.id, celltype, n.iter, q12, states )
  
  load( sprintf( "%s/%s/Figures/gg_figs_all_2.rda", ps$RESULTS_PATH, celltype ) )
  
  gg_figs.ABCDE <- ggpubr::ggarrange(
    gg_figs_all$gg.fig.5A, gg_figs_all$gg.fig.5BCD, gg_figs_all$gg.fig.5E,
    labels = c( "A", "", "E" ),
    ncol = 1, nrow = 3 )
  
  cm_to_in <- 1 / 2.54  
  ggpubr::ggexport( gg_figs.ABCDE,
                    filename = sprintf( "%s/%s_Figure_5ABCDE.pdf",
                                        path.figures, fig.filename.start ),
                    width = 21 * cm_to_in, height = 18 * cm_to_in )
  ggpubr::ggexport( gg_figs.ABCDE,
                    filename = sprintf( "%s/%s_Figure_5ABCDE.svg",
                                        path.figures, fig.filename.start ),
                    width = 21 * cm_to_in, height = 18 * cm_to_in )
}
