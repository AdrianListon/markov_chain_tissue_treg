# converting_Q_and_flow_tables.r

get.exit.flow.table <- function( ps, celltype, tissue, model.dir.name )
{
  exit.flow.file <- sprintf( 
    "%s/%s/%s/%s/parabio_fit_HDI.csv",
    ps$ANALYSIS_PATH, celltype, tissue, model.dir.name )
  if( file.exists( exit.flow.file ) ) {
    read_csv( exit.flow.file, show_col_types = FALSE ) %>% 
    dplyr::select( par, mode ) %>% 
    filter( par %in% c( "q12", "q23", "q32" ) |
      stringr::str_detect( par, "^q[1-6][1-6]_i\\[1" ) ) %>% 
    filter( !( par %in% sprintf( "q%s%s_i[1]", 1:6, 1:6 ) ) ) %>%
    mutate( par = substring( par, first = 1, last = 3 ) ) %>% 
    
    mutate( population.1 = str_match( par, "^q(\\d)[1-6]" )[ , 2 ] %>% 
              as.integer,
            population.2 = str_match( par, "^q[1-6](\\d)" )[ , 2 ] %>% 
              as.integer ) %>% 
    dplyr::rename( flow = mode ) %>% 
    mutate( flow.part = "exit" ) %>% 
    dplyr::select( par, population.1, population.2, flow, flow.part ) %>% 
    arrange( population.1, population.2 ) -> 
    exit.flow.table
    return( exit.flow.table )
  }
}


get.al_ <- function( ps, celltype, tissue, model.dir.name )
{
  al_.file <- sprintf( 
    "%s/%s/%s/%s/model_vars/al_.csv", 
    ps$ANALYSIS_PATH, celltype, tissue, model.dir.name )
  al_ <- read.csv( al_.file )$x[ 1 : 6 ]
  return( al_ )  
}


convert.exit.flow.table.into.entry.flow.table <- function( exit.flow.table, 
                                                           al_ )
{
  stopifnot( exists( "exit.flow.table" ) && !is.null( exit.flow.table ) && 
               ncol( exit.flow.table ) == 5 )
  
  exit.flow.table %>%
    mutate( al_.entry = al_[ population.1 ] ) %>%                              
    mutate( al_.exit = al_[ population.2 ] ) %>%
    mutate( al_.entry = ifelse( al_.entry == 0, 1e-9, al_.entry ) ) %>% 
    mutate( al_.exit = ifelse( al_.exit == 0, 1e-9, al_.exit ) ) %>% 
    mutate( flow = flow * ( al_.entry / al_.exit ) ) %>%
    mutate( flow.part = "entry" ) %>% 
    dplyr::rename( population.1.old = population.1,
                   population.1 = population.2 ) %>% 
    dplyr::rename( population.2 = population.1.old ) %>% 
    dplyr::select( par, population.1, population.2, flow, flow.part ) %>% 
    arrange( population.1, population.2 ) ->
    entry.flow.table
  return( entry.flow.table )
}


get.entry.flow.table <- function( ps, celltype, tissue, model.dir.name )
{
  entry.flow.table <- 
    convert.exit.flow.table.into.entry.flow.table( 
      exit.flow.table = get.exit.flow.table( 
        ps = ps, celltype = celltype, tissue = tissue, 
        model.dir.name = model.dir.name ), 
      al_ = get.al_( ps = ps, celltype = celltype, tissue = tissue,
                     model.dir.name = model.dir.name ) 
    )
  return( entry.flow.table )
}
