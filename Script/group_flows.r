# group_flows.r

read.and.clean.total.counts.data <- function( ps, celltype )
{
  read_csv( sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv",
                     ps$PROCESSED_PATH, celltype ) ) %>% 
    dplyr::rename( tissue = "Tissue",
                   total.count = "Mean",
                   total.count.sd = `Standard deviation` ) -> dTotal.counts
  return( dTotal.counts )
}


get.flows.all.tissues <- function( ps, celltype, flow.part,
                                   model.name, model.ver, 
                                   setup.mcmc.fname )
{
  TISSUES <- get_parabiosis_const()$TISSUES 
  
  if ( flow.part == "exit" ) { 
    get.flow.table <- get.exit.flow.table } else {
      get.flow.table <- get.entry.flow.table }
  tissue.s <- TISSUES[ !( TISSUES %in% c( "IEL", "LPL" ) ) ]
  model.dir <- sprintf( "%s_%i%s", model.name, 
                        mcmc.chain.n * ( mcmc.iter.n - mcmc.warmup.n ),
                        model.ver )
  flow.table.all.tissues <- NULL
  for ( i.tissue in tissue.s ) {  
    flow.table.all.tissues %>% 
      bind_rows( get.flow.table( ps = ps, celltype = celltype, tissue = i.tissue, 
                                 model.dir.name = model.dir ) %>% 
                   mutate( tissue = i.tissue ) ) -> flow.table.all.tissues }
  return( flow.table.all.tissues )
}


get.tissues.average.flows <- function( ps, celltype, flow.part,
                                       model.name, model.ver, 
                                       setup.mcmc.fname )
{
  flow.table <- get.flows.all.tissues( 
    ps = ps, celltype = celltype, flow.part = flow.part,
    model.name = model.name, model.ver = model.ver, 
    setup.mcmc.fname = setup.mcmc.fname )  
  total.counts.table <- read.and.clean.total.counts.data(                      
    ps = ps, celltype = celltype )
  
  
  flow.table %>%
    left_join( ., get_parabiosis_const()$dTissueAllOrderedGroup.f %>%                                 
                 dplyr::select( tissue.all, tissue.group, f.tissue.group ), 
               by = c( "tissue" = "tissue.all" ) ) %>%
    left_join( ., total.counts.table %>% dplyr::select( tissue, total.count ), 
               by = "tissue" ) -> tmp
  
    tmp %>% 
    # par used again technical         
    group_by( f.tissue.group, par ) %>%                                        
    mutate( sum.total.count.per.tissue.group = sum( total.count ) ) %>% 
    mutate( weight.of.tissue.in.tissue.group = 
              total.count / sum.total.count.per.tissue.group ) %>% 
    mutate( flow.mean.w.by.tissue.in.tissue.group = 
              sum( weight.of.tissue.in.tissue.group * flow ) ) %>%
    
    mutate( median.flow = median( flow ) ) %>% 
    mutate( mean.flow = mean( flow ) ) %>% 
    mutate( log10.median.flow = log10( median.flow ) ) %>% 
    mutate( log10.median.flow.1000 = log10( median.flow * 1000 ) ) %>% 
    
    arrange( par, f.tissue.group ) %>% 
    ungroup() %>% 
    dplyr::select( f.tissue.group, population.1, population.2, par, 
                   mean.flow, flow.mean.w.by.tissue.in.tissue.group, flow.part, 
                   median.flow, log10.median.flow, log10.median.flow.1000 ) %>% 
    distinct() %>% 
    
    group_by( population.1, f.tissue.group ) %>% 
    mutate( sum.flow.mean.w.per.population.1 = 
              sum( flow.mean.w.by.tissue.in.tissue.group ) ) %>% 
    ungroup() %>% 
    dplyr::select( f.tissue.group, population.1, population.2, par, 
                   mean.flow, flow.mean.w.by.tissue.in.tissue.group, flow.part, 
                   sum.flow.mean.w.per.population.1, 
                   median.flow, log10.median.flow, log10.median.flow.1000 ) ->
    average.flows.table
  return( average.flows.table )
}
