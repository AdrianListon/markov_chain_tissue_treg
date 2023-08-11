# plot_all_tissue_trajectories.r

plot_all_tissue_trajectories <- function( celltype = celltype, ps, n.iter = 0,
                                          sel_models_file_name, q12, states )
{
  pmc <- get_parabiosis_const()
  pft <- parabiosis_flow_titles( celltype = celltype )
  
  
  model.id <- gsub( pattern = ".*_([0-9]{3,4}[a-z]{0,3})_.*", replacement = "\\1",
                    x = sel_models_file_name )
  fig.filename.start <- sprintf( 
    "f%s_%s_%s_Q12=%s_%s_Mode_HDI", model.id, celltype, n.iter, q12, states )
  path.figures <- sprintf( "%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures" )
  if ( !dir.exists( path.figures ) ) { 
    dir.create( path.figures, recursive = TRUE ) }
  path.fig.filename.start <- sprintf( "%s/%s", path.figures, fig.filename.start )
  if ( !dir.exists( path.fig.filename.start ) ) {
    dir.create( path.fig.filename.start, recursive = FALSE ) }
  
  tissues <- pmc$tissue.all.ordered
  tissues <- c( tissues[ -which( tissues == "Blood" ) ], "Blood" )

  
  load_one_obj <- function( f )
  {
    env <- new.env()
    nm <- load( f, env )[ 1 ]
    return( env[[ nm ]] )
  }

  
  l_tissue_per_row <- list()
  # for( i in ( 1 : length( tissues ) ) ) {
  #   l_tissue_per_row[[ i ]] <- load_one_obj( sprintf( 
  #     "%s/f1367_Treg_20000_Q12=free_states_6_Mode_HDI_%s_trajectories_2.rda",   
  #     path.figures, tolower( tissues[ i ] ) ) )
  # }
  
  for( i in ( 1 : length( tissues ) ) ) {
    l_tissue_per_row[[ i ]] <- load_one_obj( sprintf( 
      "%s/%s_%s_trajectories_2.rda",   
      path.figures, fig.filename.start, tolower( tissues[ i ] ) ) )
  }
  

  
  n_figs_per_A4 <- 4
  if ( ( length( tissues ) %% n_figs_per_A4 ) != 0 ) {
    for ( i in ( ( length( tissues ) + 1 ) : 
                 ( ceiling( length( tissues ) / n_figs_per_A4 ) * n_figs_per_A4 ) ) ) {
      l_tissue_per_row[[ i ]] <- NULL  
    }
  }
  
  l_trajectories_multiple_pages <- list()
  for ( i in 0 : ( ceiling( length( tissues ) %/% n_figs_per_A4 ) - 1 ) ) {
    l_trajectories_multiple_pages[[ 1 + i ]] <- ggpubr::ggarrange( 
      l_tissue_per_row[[ 1 + i*n_figs_per_A4 ]], 
      l_tissue_per_row[[ 2 + i*n_figs_per_A4 ]],
      l_tissue_per_row[[ 3 + i*n_figs_per_A4 ]], 
      l_tissue_per_row[[ 4 + i*n_figs_per_A4 ]],
      ncol = 1, nrow = n_figs_per_A4, 
      labels = LETTERS[ ( 1:n_figs_per_A4 ) + i * n_figs_per_A4 ]
    )
  }
  
  l_trajectories_multiple_pages[[ 5 ]] <- ggpubr::ggarrange( 
    l_tissue_per_row[[ 17 ]], NULL, NULL, NULL,
    ncol = 1, nrow = n_figs_per_A4, 
    labels = c( LETTERS[ 17 ], "", "", "" ) )
  
  cm_to_in <- 1 / 2.54
  cm_to_in2 <- cm_to_in * 2
  ggpubr::ggexport( l_trajectories_multiple_pages,
                    filename = sprintf( "%s/Figure_S17_%s_trajectories_multipage.pdf", 
                                        path.figures, fig.filename.start ),
                    width = 21 * cm_to_in2, height = 29.7 * cm_to_in2 )

  for ( i in 0 : ( ceiling( length( tissues ) %/% n_figs_per_A4 ) ) ) {
    ggpubr::ggexport( l_trajectories_multiple_pages[[ 1 + i ]],
                      filename = sprintf( 
                        "%s/Figure S17 p%s.svg", 
                        path.figures, as.character( i + 1 ) ),
                      width = 21 * cm_to_in2, height = 29.7 * cm_to_in2 )
  }
}
