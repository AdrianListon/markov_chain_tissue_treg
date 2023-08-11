# estimate_one_markov_model.r

estimate_one_markov_model <- function( 
    celltype, ps, vector.tissues, 
    model.name, model.ver,
    use_analysis_done_as_input = FALSE, N3_modelled = TRUE,
    n.iter, sel_models_file_name, q12, states )
{                                      
  # fit markov model
  if ( use_analysis_done_as_input ) {                                            
    output.tissue.dir <- sprintf( 
      "%s/%s/%s", ps$ANALYSIS_DONE_PATH, celltype, vector.tissues[ 1 ] )
  } else {
    output.tissue.dir <- sprintf( 
      "%s/%s/%s", ps$ANALYSIS_PATH, celltype, vector.tissues[ 1 ] )     
  
    mcmc.par.set <- ( ( mcmc.iter.n - mcmc.warmup.n ) * mcmc.chain.n ) %>% 
      format( scientific = FALSE )
    output.dir <- sprintf( 
      "%s/%s_%s%s", output.tissue.dir, model.name, mcmc.par.set, model.ver ) 
    if ( !dir.exists( output.dir ) ) { 
      dir.create( output.dir, recursive = TRUE ) }
  }
  
  # LOAD_PREMODEL_VARS_AND_DATA_FROM_markov_premodel_prints
    analysis.celltype.tissue.model.vars.dir <- 
      sprintf( "%s/%s/%s/%s_%s%s/model_vars", ps$ANALYSIS_PATH, celltype, 
               vector.tissues[ 1 ], model.name, mcmc.par.set, model.ver )
    load( file = sprintf( "%s/model_vars_.RData", 
                          analysis.celltype.tissue.model.vars.dir ) )
    load( file = sprintf( "%s/pb.RData", 
                          analysis.celltype.tissue.model.vars.dir ) )
    load( file = sprintf( "%s/other_vars.RData",                                
                          analysis.celltype.tissue.model.vars.dir ) )
  
    
    read_Q_matrix <- function( output.dir.sel.chains ) {
      Q_matrix_files <- list.files( 
        output.dir.sel.chains, pattern = "^Q_matrix.*.csv", full.names = TRUE )
      stopifnot( length( Q_matrix_files ) == 1 )
      Q <- read.csv( Q_matrix_files[ 1 ] )[ , -1 ] %>% as.matrix
      rownames( Q ) <- colnames( Q )
      return( Q )
    }
    
    plot_par_densities_and_calc_Q <- function( 
    parabio.fit, sel.type, sel.chains, model.name, 
    output.dir, ps, analysis.celltype.tissue.model.vars.dir,
    iter = "all" )   # iter means number of interations AFTER warmup! 
    {                                                           
      
      output.dir <- sprintf( "%s/%s_chains_%s", output.dir, sel.type, sel.chains )
      dir.create( output.dir, recursive = TRUE )
      
      # pb.tissue, pb.hostdonor.popul, parabio.data, pb.hostdonor.tissue.popul, 
      # pb.popul
      load( file = sprintf( "%s/pb.RData",                                      
                            analysis.celltype.tissue.model.vars.dir ) )
      # other_vars = parabio.data
      load( file = sprintf( "%s/other_vars.RData",                                
                            analysis.celltype.tissue.model.vars.dir ) )
      
      
      if ( GET_FITTED_PARAMETERS_FROM_MCMC_SAMPLE_TO_R_ENV <- TRUE ) {
        sel.chains <- stringr::str_split_1( sel.chains, pattern = "_" ) %>% 
          as.numeric()
        if ( iter == "all" ) {
          # get MCMC sample    
          parabio.fit_draws <- posterior::as_draws( parabio.fit )
          draws_ok <- posterior::subset_draws( 
            parabio.fit_draws, chain = sel.chains )
          summary_draws_ok <- posterior::summarise_draws( draws_ok )
          write_csv( summary_draws_ok, 
                     file = sprintf( "%s/parabio_fit_print_nOkCh=%s.csv", 
                                     output.dir, sum( !is.na( sel.chains ) ) ) )
        }
        parabio.fit.sample <- rstan::extract( parabio.fit, permuted = FALSE )
        
        if ( GET_ALL_THE_CALCULATED_PARAMS_TO_THE_ENV <- TRUE ) {
          # get the fitted parameters ("free" parameters here only) from the MCMC 
          # sample for sel.chains
          
          if ( ASSIGN_VALUES_TO_FREE_PARAMS__BACKTICKS_R <- TRUE ) {
            for ( mp in dimnames( parabio.fit.sample )$parameters )
            {
              if ( iter == "all" ) {
                pf.sample <- as.vector( parabio.fit.sample[ , sel.chains, mp ] )
                pf.sample.mode <- mlv( pf.sample, method = "venter", type = "shorth" )
                pf.sample.density <- density( pf.sample, n = 1000, cut = 0 )
                pf.sample.density.max.x <- pf.sample.density$x[
                  which.max( pf.sample.density$y ) ]
                assign( mp, pf.sample.density.max.x )
                
                
                png( filename = sprintf( "%s/density_%s.png", output.dir, mp ),
                     width = 1024, height = 1024 )
                par( mfrow = c( 2, 1 ), mar = c( 5, 6, 1, 1 ), lwd = 3, cex.lab = 2,
                     cex.axis = 2 )
                hist( pf.sample, breaks = 30, freq = FALSE, xlab = mp, main = "" )
                abline( v = pf.sample.density.max.x, col = "blue3" )
                abline( v = pf.sample.mode, col = "red3", lty = 2 )
                plot( pf.sample.density, xlab = mp, main = "" )
                abline( v = pf.sample.density.max.x, col = "blue3" )
                abline( v = pf.sample.mode, col = "red3", lty = 2 )
                dev.off()
              } else {
                stopifnot( iter <= dim( parabio.fit.sample )[ 1 ] )
                pf.value.iter <- as.vector( parabio.fit.sample[ , sel.chains, mp ] )[ 
                  as.numeric( iter ) ]
                assign( mp, pf.value.iter )
              }
            }
          }
          
          # parse the .stan code to create the whole Q matrix in the R environment
          stan.file <- readr::read_lines( sprintf( "%s/%s.stan", ps$CODE_PATH, model.name ) )
          
          if ( FIND_STAN_CODE_AND_INITIALIZE_VECTORS__STANDARD_R <- TRUE ) {
            # get the .stan code part, which will initialize vectors in R
            # ( in .stan it is done by different code in part parameters{} )
            
            # get the .stan code part
            stan.file %>% str_which( 
              "\\/\\/ start of definitions for R", negate = FALSE ) -> n.line.definitions.start
            stan.file %>% str_which( 
              "\\/\\/ end of definitions for R", negate = FALSE ) -> n.line.definitions.end
            
            # initialize vectors
            for ( n.line in ( n.line.definitions.start + 1 ) : ( n.line.definitions.end - 1 ) ) {
              stan.file[ n.line ] %>% 
                str_replace_all( ., pattern = "\\{\\{//REMOVE IN R", 
                                 replacement = c( "{{//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "\\}\\}//REMOVE IN R", 
                                 replacement = c( "}}//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
                
                gsub( "\\/\\/\\#" , "\\#", . ) %>%                                       # remove stan comments but keep R comments
                gsub( "\\/\\/" , "", . ) %>%                                             # replace stan comments - hack to "send" some extra R code to the R parser
                
                str_replace_all( ., pattern = "\\{\\{", 
                                 replacement = c( "{{" = "XXTMP1" ) ) %>%                # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "\\}\\}", 
                                 replacement = c( "}}" = "XXTMP2" ) ) %>%                # replace .stan min ( {a, b} ) syntax
                
                str_replace_all( ., pattern = "\\{", replacement = c( "{" = "c(" ) ) %>% # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "\\}", replacement = c( "}" = ")" ) ) %>%  # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "XXTMP1", 
                                 replacement = c( "XXTMP1" = "{" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "XXTMP2", 
                                 replacement = c( "XXTMP1" = "}" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
                
                str_replace_all( ., pattern = "int", replacement = c( "int" = "" ) ) %>%
                str_replace_all( ., pattern = "real", replacement = c( "real" = "" ) ) %>%
                
                parse( text = . ) ->                                                 
                n.line.expression
              eval( n.line.expression )
            } 
          }
          
          # Assign the modes of !!VECTORED & SAMPLED!! (ONLY???) values from .stan model to R variables
          # SOURCE OF PROBLEMS - not fully automatized!
          
          if ( CONVERT_BACKTICKS_R_TO_STANDARD_R <- TRUE ) {
            
            if ( MODEL_1000p <- TRUE ) {
              for ( i in 1 : ( tn_-1 ) ) {
                if ( exists( "q14_simplex_i[1]" ) ) {                                      
                  q14_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_simplex_i[", i, "]`" ) ) %>% eval() }
              }
              
              for ( i in 1 : ( tn_-1 ) ) {
                q51_i[ i ] <- parse( text = sprintf( "%s%i%s", "`q51_i[", i, "]`" ) ) %>% eval()
                q52_i[ i ] <- parse( text = sprintf( "%s%i%s", "`q52_i[", i, "]`" ) ) %>% eval()
                q41_i[ i ] <- parse( text = sprintf( "%s%i%s", "`q41_i[", i, "]`" ) ) %>% eval()
                
                if ( exists( "q61_q65_sum_i[1]" ) ) {                                      
                  q61_q65_sum_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q61_q65_sum_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q61_q65_sum_propos_i[1]" ) ) {                                      
                  q61_q65_sum_propos_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q61_q65_sum_propos_i[", i, "]`" ) ) %>% eval() }
                
                
                
                if ( exists( "q14_min_I_i[1]" ) ) {                                      
                  q14_min_I_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_min_I_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q14_max_I_i[1]" ) ) {                                      
                  q14_max_I_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_max_I_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q14_base_i[1]" ) ) {                                      
                  q14_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_base_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q14_propos_i[1]" ) ) {                                      
                  q14_propos_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_propos_i[", i, "]`" ) ) %>% eval() }
                
                
                
                
                if ( exists( "q14_max_i[1]" ) ) {                                      
                  q14_max_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_max_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q41_u_code_i[1]" ) ) {                                      
                  q41_u_code_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q41_u_code_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q45_u_code_i[1]" ) ) {                                      
                  q45_u_code_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q45_u_code_i[", i, "]`" ) ) %>% eval() }
                
                
                if ( exists( "q36_propos_i[1]" ) ) {                                      
                  q36_propos_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_propos_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q51_u_code_i[1]" ) ) {                                      
                  q51_u_code_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q51_u_code_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q61_u_code_i[1]" ) ) {                                      
                  q61_u_code_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q61_u_code_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "al1_q14_min_i[1]" ) ) {                                      
                  al1_q14_min_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`al1_q14_min_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "al4_q41_flow_i[1]" ) ) {                                      
                  al4_q41_flow_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`al4_q41_flow_i[", i, "]`" ) ) %>% eval() }
                
                
                if ( exists( "q41_flow_i[1]" ) ) {                                      
                  q41_flow_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q41_flow_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "wgh4_free_i[1]" ) ) {                                      
                  wgh4_free_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`wgh4_free_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "wgh5_free_i[1]" ) ) {                                      
                  wgh5_free_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`wgh5_free_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "wgh6_free_i[1]" ) ) {                                      
                  wgh6_free_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`wgh6_free_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "wgh4_i[1]" ) ) {                                      
                  wgh4_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`wgh4_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "wgh5_i[1]" ) ) {                                      
                  wgh5_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`wgh5_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "wgh6_i[1]" ) ) {                                      
                  wgh6_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`wgh6_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q14_min_i[1]" ) ) {                                      
                  q14_min_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_min_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q14_i[1]" ) ) {                                      
                  q14_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q14_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_propos_i[1]" ) ) {                                      
                  q25_propos_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_propos_i[", i, "]`" ) ) %>% eval() }
                
                
                if ( exists( "q41_base_i[1]" ) ) {                                      
                  q41_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q41_base_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q41a_i[1]" ) ) {                                      
                  q41a_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q41a_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q41b_i[1]" ) ) {                                      
                  q41b_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q41b_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_sum_max_simplex_i[1]" ) ) {                                      
                  q36_sum_max_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_sum_max_simplex_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_sum_min_simplex_i[1]" ) ) {                                      
                  q36_sum_min_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_sum_min_simplex_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q36_max_S_i[1]" ) ) {                                      
                  q36_max_S_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_max_S_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_max_I_i[1]" ) ) {                                      
                  q36_max_I_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_max_I_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_max_i[1]" ) ) {                                      
                  q36_max_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_max_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_min_S_i[1]" ) ) {                                      
                  q36_min_S_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_min_S_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_min_I_i[1]" ) ) {                                      
                  q36_min_I_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_min_I_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_min_i[1]" ) ) {                                      
                  q36_min_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_min_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q36_base_i[1]" ) ) {                                              
                  q36_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_base_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_simplex_i[1]" ) ) {                                      
                  q36_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_simplex_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q36_i[1]" ) ) {                                      
                  q36_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q36_i[", i, "]`" ) ) %>% eval() }
                
                
                
                if ( exists( "q61_i[1]" ) ) {                                              
                  q61_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q61_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q63_i[1]" ) ) {                                              
                  q63_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q63_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q65_base_i[1]" ) ) {                                         
                  q65_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q65_base_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q65_i[1]" ) ) {                                              
                  q65_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q65_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q56_base_i[1]" ) ) {                                         
                  q56_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q56_base_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q56_i[1]" ) ) {                                         
                  q56_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q56_i[", i, "]`" ) ) %>% eval() }
                
                
                if ( exists( "q25_S_i[1]" ) ) {                                      
                  q25_S_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_S_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_inminmax_simplex_i[1]" ) ) {                                      
                  q25_inminmax_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_inminmax_simplex_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "al_2_q25__al_1[1]" ) ) {                                      
                  al_2_q25__al_1[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`al_2_q25__al_1[", i, "]`" ) ) %>% eval() }
                
                
                if ( exists( "q51_min_i[1]" ) ) {                                      
                  q51_min_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q51_min_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q51_max_i[1]" ) ) {                                      
                  q51_max_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q51_max_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q51_base_i[1]" ) ) {                                      
                  q51_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q51_base_i[", i, "]`" ) ) %>% eval() }
                
                
                if ( exists( "q45_max_i[1]" ) ) {                                      
                  q45_max_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q45_max_i[", i, "]`" ) ) %>% eval() }
                
                
                if ( exists( "q25_sum_max_simplex_i[1]" ) ) {                                      
                  q25_sum_max_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_sum_max_simplex_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_sum_min_simplex_i[1]" ) ) {                                      
                  q25_sum_min_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_sum_min_simplex_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q25_max_S_i[1]" ) ) {                                      
                  q25_max_S_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_max_S_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_max_I_i[1]" ) ) {                                      
                  q25_max_I_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_max_I_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_max_i[1]" ) ) {                                      
                  q25_max_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_max_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_min_S_i[1]" ) ) {                                      
                  q25_min_S_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_min_S_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_min_I_i[1]" ) ) {                                      
                  q25_min_I_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_min_I_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_min_i[1]" ) ) {                                      
                  q25_min_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_min_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q25_base_i[1]" ) ) {                                              
                  q25_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_base_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_simplex_i[1]" ) ) {                                      
                  q25_simplex_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_simplex_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q25_i[1]" ) ) {                                      
                  q25_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q25_i[", i, "]`" ) ) %>% eval() }
                
                if ( exists( "q54_base_i[1]" ) ) {                                      
                  q54_base_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q54_base_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q54_max_i[1]" ) ) {                                      
                  q54_max_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q54_max_i[", i, "]`" ) ) %>% eval() }
                if ( exists( "q54_i[1]" ) ) {                                      
                  q54_i[ i ] <- 
                    parse( text = sprintf( "%s%i%s", "`q54_i[", i, "]`" ) ) %>% eval() }
                
              }  
              if ( exists( "q12" ) ) {                                                     
                q12 <- parse( text = sprintf( "`q12`" ) ) %>% eval() 
              } else { 
                if ( q12_fixed_ == 1 ) {
                  q12 = q12_ } else { stop() }
              }
            }
          }
          
          if ( FIND_STAN_CODE_AND_ASSIGN_VALUES_TO_NONFREE_PARAMETERS__STANDARD_R <- TRUE ) {
            # get the .stan code part, which will assign values to the non-free parameters
            stan.file %>% str_which( 
              "\\/\\/ start of nonfree variables", negate = FALSE ) -> n.line.nonfree.start
            stan.file %>% str_which( 
              "\\/\\/ end of nonfree variables", negate = FALSE ) -> n.line.nonfree.end
            
            # assign the value to the non-free parameters based on the extracted .stan code
            for ( n.line in ( n.line.nonfree.start + 1 ) : ( n.line.nonfree.end - 1 ) ) {
              stan.file[ n.line ] %>%
                str_replace_all( ., pattern = "\\{\\{//REMOVE IN R", 
                                 replacement = c( "{{//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "\\}\\}//REMOVE IN R", 
                                 replacement = c( "}}//REMOVE IN R" = "" ) ) %>%         # replace .stan min ( {a, b} ) syntax
                
                gsub( "\\/\\/\\#" , "\\#", . ) %>%                                       # remove stan comments but keep R comments
                gsub( "\\/\\/" , "", . ) %>%                                             # replace stan comments - hack to "send" some extra R code to the R parser
                
                str_replace_all( ., pattern = "\\{\\{", 
                                 replacement = c( "{{" = "XXTMP1" ) ) %>%                # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "\\}\\}", 
                                 replacement = c( "}}" = "XXTMP2" ) ) %>%                # replace .stan min ( {a, b} ) syntax
                
                
                str_replace_all( ., pattern = "\\{", replacement = c( "{" = "c(" ) ) %>% # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "\\}", replacement = c( "}" = ")" ) ) %>%  # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "XXTMP1", 
                                 replacement = c( "XXTMP1" = "{" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
                str_replace_all( ., pattern = "XXTMP2", 
                                 replacement = c( "XXTMP1" = "}" ) ) %>%                 # replace .stan min ( {a, b} ) syntax
                
                str_replace_all( ., pattern = "int", replacement = c( "int" = "" ) ) %>%
                str_replace_all( ., pattern = "real", replacement = c( "real" = "" ) ) %>%
                
                parse( text = . ) ->                                                 
                n.line.expression
              eval( n.line.expression )
            } 
          }
        }
      }
      
      # having variables in environment build intensity matrix Q
      if ( HAVING_PARAMS_IN_STANDARD_R_ENV_BUILD_AND_SAVE_INTENSITY_MATRIX_Q <- TRUE ) {
        Q <- matrix( 
          0, nrow = 2 * pb.tissue.n * pb.popul.n, ncol = 2 * pb.tissue.n * pb.popul.n )
        
        stan.file %>% 
          str_which( ., "\\/\\/ start of Q alignment", negate = FALSE) -> 
          n.line.Q.alignment.start
        stan.file %>% 
          str_which( ., "\\/\\/ end of Q alignment", negate = FALSE) -> 
          n.line.Q.alignment.end
        stan.file[ ( n.line.Q.alignment.start + 1 ) : ( n.line.Q.alignment.end - 1 ) ]
        
        
        for ( n.line in ( n.line.Q.alignment.start + 1 ) : ( n.line.Q.alignment.end - 1 ) ) {
          stan.file[ n.line ] %>% 
            # replace stan comments
            gsub( "\\/\\/" , "\\#", . ) %>%
            # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\{\\{", 
                             replacement = c( "{{" = "XXTMP1" ) ) %>%                
            # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}\\}", 
                             replacement = c( "}}" = "XXTMP2" ) ) %>%                
            # replace stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\{", replacement = c( "{" = "c(" ) ) %>% 
            # replace stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "\\}", replacement = c( "}" = ")" ) ) %>%  
            # replace stan vector [a, b] syntax
            str_replace_all( ., pattern = "= \\[", 
                             replacement = c( "= \\[" = "= c(" ) ) %>%               
            str_replace_all( ., pattern = "\\];", 
                             replacement = c( "\\];" = ");" ) ) %>%
            # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "XXTMP1", 
                             replacement = c( "XXTMP1" = "{" ) ) %>%                 
            # replace .stan min ( {a, b} ) syntax
            str_replace_all( ., pattern = "XXTMP2", 
                             replacement = c( "XXTMP1" = "}" ) ) %>%                 
            parse( text = . ) ->
            n.line.expression
          eval( n.line.expression )
        } 
        
        
        dimnames( Q ) <- list( pb.hostdonor.tissue.popul, pb.hostdonor.tissue.popul )
        if ( iter == "all" ) {
          error.suffix <- ""
          if ( !( all ( Q - diag( diag( Q ) ) >= 0 ) ) ) {
            error.suffix <- paste0( error.suffix, "_Qsigns" ) }
          if ( !( all ( abs( rowSums( Q ) ) < 1e-10 ) ) ) {
            error.suffix <- paste0( error.suffix, "_rowSums0" ) }
          write.csv( Q, file = sprintf( 
            "%s/Q_matrix%s.csv", output.dir, error.suffix ), quote = FALSE )
        } else {
          write.csv( Q, file = sprintf( 
            "%s/tmp_Q_matrix%s.csv", output.dir, "_iter_x" ), quote = FALSE )
        }
      }
      return( Q )
    }  
    
    get_and_plot_trajectories <- function( 
    sel.chains, sel.type, celltype,
    output.dir, ps, analysis.celltype.tissue.model.vars.dir,
    n.iter, sel_models_file_name, q12, states )
    {
      load( file = sprintf( "%s/pb.RData",                                        
                            analysis.celltype.tissue.model.vars.dir ) )
      load( file = sprintf( "%s/other_vars.RData",                                # other_vars = parabio.data
                            analysis.celltype.tissue.model.vars.dir ) )
      output.dir.sel.chains <- sprintf( "%s/%s_chains_%s", output.dir, sel.type, sel.chains )
      
      Q <- read_Q_matrix( output.dir.sel.chains )
      
      # get limiting distribution in paired mice
      if ( FROM_Q_GET_LIMITING_DISTRIBUTION_FOR_TRAJECTORIES_PLOT <- TRUE ) {
        Qs <- Q
        n.zeroed <- 0
        Qs[ , 1 ] <- rep( 1, 2 * pb.tissue.n * pb.popul.n - n.zeroed )
        bl <- solve( t( Qs ), c( 1, rep( 0, 2 * pb.tissue.n * pb.popul.n - 1 - n.zeroed ) ) )
        names( bl )[ c( 3, pb.tissue.n * pb.popul.n + 3 ) ] <- 
          rownames( Q )[ c( 3, pb.tissue.n * pb.popul.n + 3 ) ]
      }
      
      # get average trajectories
      # INPUT (specific for sel.chains): Q
      if ( GET_AVERAGE_TRAJECTORIES <- TRUE ) {
        week.calc <- seq( w0i_, pb.figure.week.max, length.out = 200 )
        alpha.calc <- t( sapply( week.calc, function( w )
          as.vector( a0i_ %*% expm( 7 * ( w - w0i_ ) * Q ) ) ) )
        y.calc <- lapply( pb.tissue, function( tis ) {
          tissue.idx <- grep( sprintf( "\\.%s\\.", tis ), pb.hostdonor.tissue.popul )
          yc <- alpha.calc[ , tissue.idx ]
          yc <- sweep( yc, 1, rowSums( yc ), "/" )
          colnames( yc ) <- pb.hostdonor.popul
          yc
        } )
        names( y.calc ) <- pb.tissue
      }
      
      # plot trajectories
      if ( PLOT_TRAJECTORIES <- TRUE ) {
        l_ggfigs <- list()
        
        for ( tis in pb.tissue )
        { 
          l_ggfigs[[ which( pb.tissue == tis ) ]] <- list()
          names( parabio.data ) <- gsub( 
            pattern = "celltype", replacement = tolower( celltype ), 
            names( parabio.data ) )
          tissue.week <- parabio.data[ parabio.data$tissue == tis, "week" ]
          tissue.data <- parabio.data[ parabio.data$tissue == tis,
                                       pb.hostdonor.popul ]
          tissue.data <- sweep( tissue.data, 1, rowSums( tissue.data ), "/" )
          
          limit.data <- bl[ grep( sprintf( "\\.%s\\.", tis ),
                                  pb.hostdonor.tissue.popul ) ]
          limit.data <- limit.data / sum( limit.data )
          names( limit.data ) <- pb.hostdonor.popul
          
          for ( pop in pb.popul )
          {
            host.pop <- sprintf( "host.%s", pop )
            donor.pop <- sprintf( "donor.%s", pop )
            ggdata.calc <- data.frame( x = week.calc,
                                       y1 = y.calc[[ tis ]][ , host.pop ],
                                       y2 = y.calc[[ tis ]][ , donor.pop ] )
            ggdata.exp <- data.frame( x = tissue.week,
                                      y1 = tissue.data[[ host.pop ]], 
                                      y2 = tissue.data[[ donor.pop ]] )
            ggdata.exp.mean <- data.frame( t( sapply( pb.week, function( w )
              c( w, mean( tissue.data[ tissue.week == w, host.pop ] ),
                 mean( tissue.data[ tissue.week == w, donor.pop ] ) ) ) ) )
            names( ggdata.exp.mean ) <- c( "x", "y1", "y2" )
            y.max <- 1
            
            ggfig <- ggplot( ggdata.calc ) +
              geom_hline( yintercept = limit.data[ host.pop ], linetype = "dashed" ) +
              geom_hline( yintercept = limit.data[ donor.pop ], linetype = "dotted" ) +
              geom_line( aes( x, y1 ), color = "blue3" ) +
              geom_line( aes( x, y2 ), color = "red3" ) +
              geom_point( aes( x, y1 ), data = ggdata.exp, size = 1,
                          position = position_nudge( - 0.1 ), color = "blue3" ) +
              geom_point( aes( x, y2 ), data = ggdata.exp, size = 1,
                          position = position_nudge( 0.1 ), color = "red3" ) +
              geom_point( aes( x, y1 ), data = ggdata.exp.mean, size = 3,
                          position = position_nudge( - 0.1 ), color = "blue3" ) +
              geom_point( aes( x, y2 ), data = ggdata.exp.mean, size = 3,
                          position = position_nudge( 0.1 ), color = "red3" ) +
              scale_x_continuous( breaks = seq( 0, pb.figure.week.max, by = 2 ),
                                  name = "Week" ) +
              scale_y_continuous( 
                limits = c( 0, y.max ), 
                name = sprintf( "Fraction in %s %s", pb.tissue.label[ tis ], 
                                pb.parent.popul.label ) ) +
              ggtitle( sprintf( 
                "%s %s", pb.tissue.label[ tis ],
                # renaming for Tissue Treg paper
                ifelse( ( pb.popul.label[ pop ] == "Treg Naive" ), 
                        "Treg Resting", pb.popul.label[ pop ] ) ) ) +           
              theme_bw() +
              theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text( size = 14 ),
                     axis.title = element_text( size = 14 ),
                     plot.title = element_text( size = 16 ) )
            
            ggsave( sprintf( "%s/%s_%s.png", output.dir.sel.chains, tis, pop ), ggfig, 
                    width = 8, height = 6 )
            l_ggfigs[[ which( pb.tissue == tis ) ]][[ which( pb.popul == pop ) ]] <- ggfig
          }
        }
        
        # all_trajectories figure
        if ( length( pb.tissue ) == 4 ) {
          ggfigs_ggpubr <- ggpubr::ggarrange( 
            # naive
            l_ggfigs[[ 1 ]][[ 1 ]], l_ggfigs[[ 2 ]][[ 1 ]], l_ggfigs[[ 3 ]][[ 1 ]], l_ggfigs[[ 4 ]][[ 1 ]],
            # activated
            l_ggfigs[[ 1 ]][[ 2 ]], l_ggfigs[[ 2 ]][[ 2 ]], l_ggfigs[[ 3 ]][[ 2 ]], l_ggfigs[[ 4 ]][[ 2 ]],
            # cd69+
            l_ggfigs[[ 1 ]][[ 3 ]], l_ggfigs[[ 2 ]][[ 3 ]], l_ggfigs[[ 3 ]][[ 3 ]], l_ggfigs[[ 4 ]][[ 3 ]],
            labels = LETTERS[ 1 : ( 3 * length( pb.tissue ) ) ],
            ncol = length( pb.tissue ), nrow = 3 )  
        } else {
          ggfigs_ggpubr <- ggpubr::ggarrange( 
            # naive
            l_ggfigs[[ 1 ]][[ 1 ]], l_ggfigs[[ 2 ]][[ 1 ]], l_ggfigs[[ 3 ]][[ 1 ]], 
            # activated
            l_ggfigs[[ 1 ]][[ 2 ]], l_ggfigs[[ 2 ]][[ 2 ]], l_ggfigs[[ 3 ]][[ 2 ]], 
            # cd69+
            l_ggfigs[[ 1 ]][[ 3 ]], l_ggfigs[[ 2 ]][[ 3 ]], l_ggfigs[[ 3 ]][[ 3 ]], 
            labels = LETTERS[ 1 : ( 3 * length( pb.tissue ) ) ],
            ncol = length( pb.tissue ), nrow = 3 )  
        }
        
        
        
        
        ggsave( sprintf( "%s/%s.png", output.dir.sel.chains, "all_trajectories" ), ggfigs_ggpubr,
                width = 29.7*2.5, height = 21.0*2.5, units = "cm" )
        
        model.id <- gsub( pattern = ".*_([0-9]{3,4}[a-z]{0,3})_.*", replacement = "\\1",
                          x = sel_models_file_name )
        fig.filename.start <- sprintf( 
          "f%s_%s_%s_Q12=%s_%s_Mode_HDI", model.id, celltype, n.iter, q12, states )
        path.figures <- sprintf( "%s/%s/%s", ps$RESULTS_PATH, celltype, "Figures" )
        
        
        ggsave( sprintf( "%s/%s_all_trajectories_%s.png",
                         path.figures, fig.filename.start, pb.tissue[ 2 ] ),
                ggfigs_ggpubr, width = 29.7*2.5, height = 21.0*2.5, units = "cm" )
        
        
        
        gg_tissue_i_trajectories <- ggpubr::ggarrange( 
          l_ggfigs[[ 2 ]][[ 1 ]], l_ggfigs[[ 2 ]][[ 2 ]], l_ggfigs[[ 2 ]][[ 3 ]],
          ncol = 3, nrow = 1 )
        gg_tissue_i_trajectories %>% save( ., file = sprintf( 
          "%s/%s_%s_trajectories_2.rda",
          path.figures, fig.filename.start, pb.tissue[ 2 ] ) )
        
        if ( pb.tissue[ 2 ] == "brain" ) {
          gg_tissue_blood_trajectories <- ggpubr::ggarrange( 
            l_ggfigs[[ 1 ]][[ 1 ]], l_ggfigs[[ 1 ]][[ 2 ]], l_ggfigs[[ 1 ]][[ 3 ]],
            ncol = 3, nrow = 1 )
          gg_tissue_blood_trajectories %>% save( ., file = sprintf( 
            "%s/%s_%s_trajectories_2.rda",
            path.figures, fig.filename.start, "blood" ) )
        }
      }
    }  
    
    test_eqeq <- function(
    sel.chains, sel.type, celltype, 
    output.dir, ps, analysis.celltype.tissue.model.vars.dir, iter = NA,
    write_csv = TRUE )
    {
      tol <- 1e-7
      output.dir.sel.chains <- sprintf( 
        "%s/%s_chains_%s", output.dir, sel.type, sel.chains )
      
      if ( !is.na( iter ) ) {
        Q <- plot_par_densities_and_calc_Q(
          parabio.fit = parabio.fit, sel.type = sel.type, sel.chains = sel.chains, 
          model.name = model.name, output.dir = output.dir, ps = ps, 
          analysis.celltype.tissue.model.vars.dir = 
            analysis.celltype.tissue.model.vars.dir, iter = iter ) 
      } else {
        
        Q <- read_Q_matrix( output.dir.sel.chains )  
      }    
      
      load( file = sprintf( "%s/other_vars.RData",                                
                            analysis.celltype.tissue.model.vars.dir ) )
      
      eqeq_table <- tibble( 
        node = 1 : dim( Q ), influx = NA, outflux = NA, equal = NA, iter = iter )
      
      Qd0 <- Q - diag( diag( Q ) )
      for ( i in ( 1 : dim( Q )[ 1 ] ) ) {
        eqeq_table$outflux[ i ] <- -c( al_, al_ )[ i ] * Q[ i, i ]
        eqeq_table$influx[ i ] <- t( c( al_, al_ ) ) %*% Qd0[ , i ]
        eqeq_table$equal[ i ] <- abs( 
          eqeq_table$outflux[ i ] - eqeq_table$influx[ i ] ) < tol
      }
      if ( write_csv ) { write.csv( eqeq_table, sprintf( 
        "%s/aa_eqeq_table_iter_%s.csv", output.dir.sel.chains,
        ifelse( is.na( iter), "all", iter ) ) ) }
      return( all( eqeq_table$equal ) )
    }
    
    get_Qcounts <- function(
    sel.chains, sel.type, celltype,
    output.dir, ps, analysis.celltype.tissue.model.vars.dir ) {
      load( file = sprintf( "%s/pb.RData",                                        
                            analysis.celltype.tissue.model.vars.dir ) )
      # other_vars = parabio.data
      load( file = sprintf( "%s/other_vars.RData",                                
                            analysis.celltype.tissue.model.vars.dir ) )
      output.dir.sel.chains <- sprintf( 
        "%s/%s_chains_%s", output.dir, sel.type, sel.chains )
      
      read_csv( sprintf( "%s/Total_counts/parabiosis_model_input_%s_counts.csv",
                         ps$PROCESSED_PATH, celltype ) ) -> d.celltype.counts
      f_al_to_ncells <- ( d.celltype.counts[ 
        d.celltype.counts$Tissue == "Blood", "Mean" ] %>% 
          unlist %>% as.numeric ) / 
        sum( al_[ 1 : 3 ] )
      
      Q <- read_Q_matrix( output.dir.sel.chains )
      
      Qcounts <- Q - diag( diag( Q ) )
      Qcounts_abs <- Qcounts
      for ( i in ( 1 : dim( Qcounts )[ 1 ] ) ) { 
        Qcounts_abs[ i, ] <- c( al_, al_ )[ i ] * Qcounts[ i, ] * f_al_to_ncells
        Qcounts[ i, ] <- c( al_, al_ )[ i ] * Qcounts[ i, ]
      }  
      diag( Qcounts ) <- NA
      diag( Qcounts_abs ) <- NA
      write.csv( Qcounts, file = sprintf( 
        "%s/Qcounts_matrix.csv", output.dir.sel.chains ), quote = FALSE )
      write.csv( Qcounts_abs, file = sprintf(
        "%s/Qcounts_abs_matrix.csv", output.dir.sel.chains ), quote = FALSE )
      
      return( Qcounts )
    }
    

  # decide whether to run new MCMC simulations  
  MCMC_is_done <- file.exists( sprintf( 
    "%s/parabio_fit%i.rda", output.dir, mcmc.iter.n ) )                         
  if ( MCMC_is_done ) {
    load( sprintf( "%s/parabio_fit%i.rda", output.dir, mcmc.iter.n ) )
  } else {
    # run new MCMC simulations
    
    # copy model first to avoid access from multiple sites ( e.g., on cluster )
    file.copy( from = sprintf( "%s/%s.stan", ps$CODE_PATH, model.name ),
               to = sprintf( "%s/%s.stan", output.dir, model.name ) )
    cat( "\nTranslating from Stan to C++ and compiling in C++ ...\n" )
    
    # translate from Stan to C++ and compile in C++
    parabio.model <- stan_model( 
      sprintf( "%s/%s.stan", output.dir, model.name ), verbose = FALSE )        
    
    # estimate the model by new MCMC simulation:
    cat( "\n############\n" )
    cat( "New MCMC simulation has started.\n" )
    cat( "############\n\n" )
    data_list <- list( host_donor_rate_ = host_donor_rate_,
                       use_wghs_ = use_wghs_, 
                       use_hdr2_ = use_hdr2_, hdr2_lower_ = hdr2_lower_,
                       max_flow_to_N3_ = max_flow_to_N3_,
                       q12_ = q12_, q12_fixed_ = q12_fixed_, 
                       b_rate_rel_to_blood_naive_ = b_rate_rel_to_blood_naive_,                   
                       n_ = n_, wn_ = wn_, pn_ = pn_, tn_ = tn_,
                       w0i_ = w0i_, w_ = w_, wi_ = wi_, ti_ = ti_, x_ = x_,    
                       al_ = al_, a0i_ = a0i_ )
    model.sim <- purrr::quietly( sampling )(                                   
      diagnostic_file = sprintf( "%s/diagnostic_file.csv", output.dir ),
      sample_file = sprintf( "%s/sample_file.csv", output.dir ),
      verbose = TRUE,
      object = parabio.model, data = data_list, 
      iter = mcmc.iter.n, warmup = mcmc.warmup.n, chains = mcmc.chain.n,        
      control = list( adapt_delta = sampling.adapt_delta,
                      max_treedepth = 10 ) )                                    
    # MCMC simulation finished
    
    parabio.fit <- model.sim$result
    save( parabio.fit, file = sprintf( 
      "%s/parabio_fit%i.rda", output.dir, mcmc.iter.n ) )                       
    # MCMC simulation saved
    
    # save warnings, output, and messages
    sink( sprintf( "%s/parabio_fit_sampling_and_warnings.txt", output.dir ) )
    cat( "Warning messages:\n", model.sim$warnings, "\n" )
    cat( "Output:\n",           model.sim$output,   "\n" )
    cat( "Messages:\n",         model.sim$messages, "\n" )
    sink()
  }
  
  if ( PRINT_PARABIO_FIT_ALL_CHAINS <- TRUE ) {
    sink( file = sprintf( "%s/parabio_fit_print.txt", output.dir ) )
    width.backup <- getOption( "width" ); options( width = 200 ) 
    print( parabio.fit )
    options( width = width.backup )
    sink()
  }
  
  
  # Which chains to get results for:
  sel.chains_set_max <- select_mcmc_chains(
    parabio.fit = parabio.fit, mcmc.chain.n = mcmc.chain.n, type = "max" ) %>%
    paste( collapse = "_" )
  sel.chains_set_3mad <- select_mcmc_chains(
    parabio.fit = parabio.fit, mcmc.chain.n = mcmc.chain.n, type = "3mad" ) %>%
    paste( collapse = "_" )
  sel.chains.table <- tibble( 
    sel.chains_set = c( 1 : mcmc.chain.n, sel.chains_set_max, 
                        sel.chains_set_3mad ),
    sel.chains_type = c( rep( "each", mcmc.chain.n ), "max", "3mad" ) )
  
  write.csv( sel.chains.table, sprintf( 
    "%s/sel_chains_table.csv", output.dir ) )
  
  # either all, which takes ~5 minutes:
  # for( i in 1 : nrow( sel.chains.table ) )  {

  # or only max to save time:
  for( i in ( which( sel.chains.table$sel.chains_type == "max" ) ) ) {
      plot_par_densities_and_calc_Q(
      sel.chains = sel.chains.table$sel.chains_set[ i ],
      sel.type = sel.chains.table$sel.chains_type[ i ],
      parabio.fit = parabio.fit,
      model.name = model.name, 
      output.dir = output.dir, ps = ps, 
      analysis.celltype.tissue.model.vars.dir = 
        analysis.celltype.tissue.model.vars.dir ) ->                   
      Q  
    
    test_eqeq( 
      sel.chains = sel.chains.table$sel.chains_set[ i ],
      sel.type = sel.chains.table$sel.chains_type[ i ],
      celltype = celltype,
      output.dir = output.dir, ps = ps,
      analysis.celltype.tissue.model.vars.dir =
        analysis.celltype.tissue.model.vars.dir ) ->
      final_fit_eqeq
    
    # if debugging needed write where how the equilibrium is not met
    if ( !final_fit_eqeq & sel.chains.table$sel.chains_set[ i ] == "1" & 
         sel.chains.table$sel.chains_type[ i ] == "each" ) {
      for ( iter in sample( 1 : 100, 10, replace = FALSE ) )
      {
        test_eqeq( 
        sel.chains = sel.chains.table$sel.chains_set[ i ],
        sel.type = sel.chains.table$sel.chains_type[ i ],
        celltype = celltype, output.dir = output.dir, ps = ps, 
        analysis.celltype.tissue.model.vars.dir =
          analysis.celltype.tissue.model.vars.dir, iter = iter )
      }
    }
      
    get_Qcounts( 
      sel.chains = sel.chains.table$sel.chains_set[ i ],
      sel.type = sel.chains.table$sel.chains_type[ i ],
      celltype = celltype, output.dir = output.dir, ps = ps, 
      analysis.celltype.tissue.model.vars.dir =
        analysis.celltype.tissue.model.vars.dir )
    
    get_and_plot_trajectories(
      sel.chains = sel.chains.table$sel.chains_set[ i ],
      sel.type = sel.chains.table$sel.chains_type[ i ],
      celltype = celltype, output.dir = output.dir, ps = ps, 
      analysis.celltype.tissue.model.vars.dir =
        analysis.celltype.tissue.model.vars.dir,
      n.iter = n.iter, sel_models_file_name = sel_models_file_name, 
      q12 = q12, states = states )
  }
}
