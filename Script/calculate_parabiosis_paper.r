# calculate_parabiosis_paper.r
 
# Copyright (c) 2023 University of Cambridge and Babraham Institute (United Kingdom)
#
# Software written by Vaclav Gergelits and Carlos P. Roca as research funded by the European Union and the United Kingdom
#
# This software may be modified and distributed under the terms of the MIT
# license. See the LICENSE file for details.

# Runs MCMC simulation of parabiosis study of T cell kinetics, creating all 
# figures and tables used in tissue Treg manuscript.


calculate_parabiosis_paper <- function( 
    CELLTYPE = "Treg", 
    F_TISSUES = NA,
    VECTOR_TISSUES_0 = c( "Spleen" ),
    MODEL_NAME = "markov_model__0001",
    Q12_IS_FIXED = 0, 
    Q12_VALUE = -0.5,
    F_MODEL_VER_EXT = NA,
    STATES_5 = FALSE,
    setup.mcmc.fname = "setup_simulation_parameters.r",
    F_PS = PS,
    changed_paramaters_csv = NULL,
    use_analysis_done_as_input = FALSE,
    PB_WEEK0I = 0, 
    EPS = 0.00595,
    HOST_DONOR_RATE_ = NA,
    NORMALIZE = TRUE,
    USE_WGHS_ = 1, 
    USE_HDR2_ = 1, 
    BRF = 1.1,
    MAX_FLOW_TO_N3_ = 10 )
{
  if( VECTOR_TISSUES_0 == "NULL" ) VECTOR_TISSUES_0 <- NULL
  # FOR REMOTE SESSIONS OR WINDOWS (when not opened within RStudio Project):
  # if ( ! is.null( "F_PS" ) ) { PS = F_PS } else {
  #   if( grepl( pattern = "^C\\:\\/", x = getwd() ) ) { source( "SET_PATHS.r" ) } else
  #   { source( "Dropbox (Cambridge University)/BI_20220810/d/AAP/PMM_GitHub/SET_PATHS.r" ) }
  # }
  
  suppressPackageStartupMessages( {
    library( digest )
    library( expm )
    library( gridExtra )   
    library( gtools )
    library( modeest )
    library( rstan )
    library( stringr )     
    library( tidyverse )
  } )
  # can be used depending on the tidyverse / readr library version loaded
  try( options( readr.show_col_types = FALSE ) )  
  
  # the package posterior needed but better to avoid loading its functions 
  # and masking the default ones ( e.g., mad() )
  # similar for ggpubr
  stopifnot( nchar( system.file( package = "posterior" ) ) > 0 ) 
  stopifnot( nchar( system.file( package = "ggpubr" ) ) > 0 )
  
  # load my own functions
  source( sprintf( "%s/%s", PS$CODE_PATH, "load_my_functions.r" ) )
  load_my_functions( ps = PS )
  
  # create directories
  for ( i in 1 : length( PS ) ) {
    if ( !dir.exists( PS[[ i ]] ) ) { 
      dir.create( PS[[ i ]], recursive = TRUE ) } 
    }
  
  
  # SETUP PARAMETERS
  pmc <- get_parabiosis_const()                                                 
  if ( !is.na( F_TISSUES ) ) { TISSUES = F_TISSUES } else { 
    TISSUES = pmc$TISSUES } 
  N_TISSUES = length( VECTOR_TISSUES_0 ) + 1;           
  if ( Q12_IS_FIXED != 0 ) { MODEL_NAME = sprintf( "%s_q12F", MODEL_NAME ) }
  MODEL_VER <- sprintf( 
    "_tn_=%s__Q12=%s", N_TISSUES,                           
    ifelse( Q12_IS_FIXED == 1, ceiling( Q12_VALUE * 1e4 ) / 1e4, "free" ) ); 
  if ( !is.na( F_MODEL_VER_EXT ) ) { 
    MODEL_VER = sprintf( "%s%s", MODEL_VER, F_MODEL_VER_EXT ) }

  
  source( sprintf( "%s/%s", PS$CODE_PATH, setup.mcmc.fname ) )
  SEL_MODELS_FILE_NAME <- sprintf( 
    "%s_selected_models_vX_%s_%s_Q12=%s_%s.csv", CELLTYPE, 
    mcmc.chain.n * ( mcmc.iter.n - mcmc.warmup.n ),
    substring( MODEL_NAME, first = nchar( "markov_model__" ) + 1 ), 
    ifelse( Q12_IS_FIXED == 1, ceiling( Q12_VALUE * 1e4 ) / 1e4, "free" ),
    sprintf( "states_%s", ifelse( STATES_5, "5", "6" ) ) )

  
  
  # Compute new MCMC simulation  
  for ( TIS_I in TISSUES ) { 
    VECTOR_TISSUES <- c( TIS_I, VECTOR_TISSUES_0 )
    if( length( VECTOR_TISSUES ) > length( unique( VECTOR_TISSUES ) ) ) {                         
      VECTOR_TISSUES <- c( unique( VECTOR_TISSUES ), "MLN" ) }
    COMPLEM_TISSUES = pmc$TISSUES[ !( pmc$TISSUES %in% VECTOR_TISSUES ) ]    
    
    markov_premodel_prints_eps( 
      celltype = CELLTYPE, ps = PS, vector.tissues = VECTOR_TISSUES, 
      q12_fixed_ = Q12_IS_FIXED, q12_ = Q12_VALUE, 
      model.name = MODEL_NAME, model.ver = MODEL_VER, states_5 = STATES_5,
      complem_tissues = COMPLEM_TISSUES,
      f_pb.week0i = PB_WEEK0I, f_eps = EPS,
      f_host_donor_rate_ = HOST_DONOR_RATE_,
      f_use_wghs_ = USE_WGHS_, f_use_hdr2_ = USE_HDR2_, f_brf = BRF,
      f_max_flow_to_N3_ = MAX_FLOW_TO_N3_,
      NORMALIZE = NORMALIZE ) 
    
    estimate_one_markov_model( 
      celltype = CELLTYPE, ps = PS, vector.tissues = VECTOR_TISSUES,
      model.name = MODEL_NAME, model.ver = MODEL_VER,
      use_analysis_done_as_input = use_analysis_done_as_input,
      n.iter = mcmc.chain.n * ( mcmc.iter.n - mcmc.warmup.n ),
      sel_models_file_name = SEL_MODELS_FILE_NAME,
      q12 = ifelse( Q12_IS_FIXED == 1, ceiling( Q12_VALUE * 1e4 ) / 1e4, "free" ),
      states = sprintf( "states_%s", ifelse( STATES_5, "5", "6" ) ) )          
    
    while ( sink.number() > 0 ) sink()                                          
  }
  
  
  # Calculate the MCMC accuracy:
  for ( TIS_I in pmc$TISSUES ) {
    calculate_HDI( 
      tissue = TIS_I, celltype = CELLTYPE, ps = PS,             
      model.name = MODEL_NAME, model.ver = MODEL_VER,
      setup.mcmc.fname = setup.mcmc.fname,
      use_analysis_done_as_input = use_analysis_done_as_input, 
      sel_models_file_name = SEL_MODELS_FILE_NAME, 
      create_new_csv = !is.na( SEL_MODELS_FILE_NAME ) & 
        ( TIS_I == pmc$TISSUES[ 1 ] ) )                         
  }
  
  visualize_results( 
    celltype = CELLTYPE, ps = PS,
    n.iter = mcmc.chain.n * ( mcmc.iter.n - mcmc.warmup.n ),
    sel_models_file_name = SEL_MODELS_FILE_NAME,
    q12 = ifelse( Q12_IS_FIXED == 1, ceiling( Q12_VALUE * 1e4 ) / 1e4, "free" ),  
    states = sprintf( "states_%s", ifelse( STATES_5, "5", "6" ) ),
    model.name = MODEL_NAME, model.ver = MODEL_VER, 
    setup.mcmc.fname = setup.mcmc.fname )  
  
  sessionInfo()
}
