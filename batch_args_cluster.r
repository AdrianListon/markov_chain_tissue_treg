# batch_args_cluster.r

# set the paths of the project
FILL_IN_MY_PROJECT_PATH <- ""
MY_PROJECT_PATH <- ifelse( nchar( FILL_IN_MY_PROJECT_PATH ) != 0,
                           FILL_IN_PROJECT_PATH, getwd() )
source( sprintf( "%s/SET_PATHS.r", MY_PROJECT_PATH ) )


# set cell type and tissue to estimate the Markov Model for
FF_CELLTYPE <- "Treg"
FF_TISSUE <- "Pancreas"
FF_TIS_0 <- "NULL"
args <- commandArgs( TRUE )
if ( length( args ) >= 1 ) {
  FF_CELLTYPE <- args[[ 1 ]]
  if ( length( args ) >= 2 ) {
    FF_TISSUE <- args[[ 2 ]]
    if ( length( args ) == 3 ) {
      FF_TIS_0 <- args[[ 3 ]]
    }
  }
}

source( sprintf( "%s/%s", CODE_PATH, "calculate_parabiosis_paper.r" ) )
source( sprintf( "%s/%s", CODE_PATH, "set_model_parameters.r" ) )

calculate_parabiosis_paper( 
  CELLTYPE = FF_CELLTYPE, 
  F_TISSUES = FF_TISSUE,
  VECTOR_TISSUES_0 = FF_TIS_0,       
  MODEL_NAME = "markov_model__0001",
  Q12_IS_FIXED = 0, 
  Q12_VALUE = Q12_VALUE, 
  F_MODEL_VER_EXT = "v01", 
  setup.mcmc.fname = "setup_simulation_parameters.r",
  PB_WEEK0I = WEEK0, 
  HOST_DONOR_RATE_ = HDR, 
  USE_WGHS_ = U_WGHS, 
  USE_HDR2_ = U_HDR2, 
  BRF = BRF, 
  MAX_FLOW_TO_N3_ = MAX_FLOW_TO_N3_FACTOR )
