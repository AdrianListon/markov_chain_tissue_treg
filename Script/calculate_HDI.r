# calculate_HDI.r

calculate_HDI <- function( tissue, celltype, ps,
                           model.name, model.ver,
                           setup.mcmc.fname,
                           use_analysis_done_as_input = FALSE,
                           sel_models_file_name = NA, create_new_csv )
{
  if ( use_analysis_done_as_input ) {
    output.tissue.dir <- 
      sprintf( "%s/%s/%s", ps$ANALYSIS_DONE_PATH, celltype, tissue )
  } else { output.tissue.dir <- 
      sprintf( "%s/%s/%s", ps$ANALYSIS_PATH, celltype, tissue ) }
  
  model.dir <- sprintf( "%s_%i%s", model.name, 
                        mcmc.chain.n * ( mcmc.iter.n - mcmc.warmup.n ),
                        model.ver )
  output.dir <- sprintf( "%s/%s", output.tissue.dir, model.dir )
  
  if ( file.exists( sprintf( "%s/parabio_fit_HDI.csv", output.dir ) ) ) {
    cat( sprintf( 
      "\nparabio_fit_HDI.csv already exists for %s.\nNo HDI calculation now.\n", 
      tissue ) )
  } else {
    
  if ( file.exists( sprintf( "%s/parabio_fit%i.rda", output.dir, mcmc.iter.n ) ) & 
       ( file.exists( sprintf( "%s/parabio_fit_print.txt", output.dir ) ) ) & 
       ( file.info( sprintf( "%s/parabio_fit_print.txt", output.dir ) )$size > 2000 ) ) {
    print( sprintf( 
      "%s -- parabio_fitXXXX.rda file exists and pararbio_fit_print.txt ok", tissue ) )
    load( sprintf( "%s/parabio_fit%i.rda", output.dir, mcmc.iter.n ) )
    
    parabio.fit.chain <- select_mcmc_chains( 
      parabio.fit = parabio.fit, mcmc.chain.n = mcmc.chain.n, type = "max" )
    
    parabio.fit.sample <- rstan::extract( parabio.fit, permuted = FALSE )
    
    dParamsHDI <- NULL
    for ( mp in dimnames( parabio.fit.sample )$parameters )
    {
      pf.sample <- as.vector( parabio.fit.sample[ , parabio.fit.chain, mp ] )
      pf.sample.mode <- mlv( pf.sample, method = "venter", type = "shorth" )
      
      HDI.80.50 <- get_HDI_4_values( 
        pf.sample, credMass.big = 0.80, credMass.small = 0.5 )
      
      dParamsHDI %>% 
        bind_rows( tibble( 
          par = mp, mode = pf.sample.mode, 
          hdi.80.lower = HDI.80.50[ 1 ] %>% as.numeric, 
          hdi.80.upper = HDI.80.50[ 4 ] %>% as.numeric,
          hdi.50.lower = HDI.80.50[ 2 ] %>% as.numeric, 
          hdi.50.upper = HDI.80.50[ 3 ] %>% as.numeric,
          mean = mean( pf.sample, na.rm = TRUE ),
          q50 =   quantile( pf.sample, 0.500 ),
          q2.5 =  quantile( pf.sample, 0.025 ),
          q97.5 = quantile( pf.sample, 0.975 ),
          q10 =   quantile( pf.sample, 0.100 ),
          q90 =   quantile( pf.sample, 0.900 ),
          q25 =   quantile( pf.sample, 0.250 ),
          q75 =   quantile( pf.sample, 0.750 ) ) ) -> 
        dParamsHDI
    }
    dParamsHDI %>% mutate( 
      hdi.80.ratio = hdi.80.upper / hdi.80.lower,
      hdi.50.ratio = hdi.50.upper / hdi.50.lower ) ->
      dParamsHDI
    dParamsHDI %>% write_csv( ., sprintf( "%s/parabio_fit_HDI.csv", output.dir ) )
    
    if ( !is.na( sel_models_file_name ) ) {
      dSelected_models <- tibble( 
        celltype = celltype, tissue = tissue, 
        model.id = str_match( model.dir, "markov_model__(.*?)_")[ , 2 ],
        model.dir.name = model.dir )
    } 
  
    } else { 
    if ( !is.na( sel_models_file_name ) ) {
      dSelected_models <- tibble( 
        celltype = celltype, tissue = tissue, 
        model.id = NA, model.dir.name = NA )
    }
  }
  
  sel_models_file_name_full <- 
    sprintf( "%s/%s", ps$ANALYSIS_PATH, sel_models_file_name )
  dSelected_models %>% 
    write.table( sel_models_file_name_full, sep = ",", row.names = FALSE,
                 col.names = !file.exists( 
                   sel_models_file_name_full ) | create_new_csv, 
                 append = !create_new_csv )
  }
}
