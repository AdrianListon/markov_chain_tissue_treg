# preproc_fs_for_figures_2_4_3_5.r

get.dAllTissues_Models_HDI_simple <- function( celltype, selected.models.file )
{
  selected.models.table <- read_csv( selected.models.file )
  
  dAllTissues_Models_HDI <- NULL
  for ( i in 1 : nrow( selected.models.table ) ) {
    i.tissue <- selected.models.table$tissue[ i ]
    i.model.dir.name <- selected.models.table$model.dir.name[ i ]
    if ( !is.na( i.model.dir.name ) ) {
      dAllTissues_Models_HDI %>%

        bind_rows( read_csv( sprintf( "%s/%s/%s/%s/parabio_fit_HDI.csv",
                                      ANALYSIS_PATH, celltype, 
                                      i.tissue, i.model.dir.name ) ) %>%
                     mutate( celltype = celltype,
                             tissue = i.tissue,
                             model.dir.name = i.model.dir.name ) %>%
                     dplyr::select( celltype, tissue, model.dir.name, 
                                    everything() ) ) ->
        dAllTissues_Models_HDI
    } else { print( sprintf( "WARNING: The tissue %s is missing", i.tissue ) ) }
  }
  return( dAllTissues_Models_HDI )
}

