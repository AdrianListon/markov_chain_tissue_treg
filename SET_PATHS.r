# SET_PATHS.r

PROJECT_PATH = MY_PROJECT_PATH;

CODE_PATH = sprintf( '%s/%s', PROJECT_PATH, 'Script' );
RAW_PATH = sprintf( '%s/%s', PROJECT_PATH, 'Raw' );                                                 
PROCESSED_PATH = sprintf( '%s/%s', PROJECT_PATH, 'Processed' );             
ANALYSIS_PATH = sprintf( '%s/%s', PROJECT_PATH, 'Analyzed' );
ANALYSIS_DONE_PATH = sprintf( '%s/%s', PROJECT_PATH, 'Analyzed_ZENODO' );
RESULTS_PATH = sprintf( '%s/%s', PROJECT_PATH, 'Results' );
PS <- list( PROJECT_PATH, CODE_PATH, RAW_PATH, 
            PROCESSED_PATH, ANALYSIS_PATH, RESULTS_PATH,
            ANALYSIS_DONE_PATH )
names( PS ) <- c( 'PROJECT_PATH', 'CODE_PATH', 'RAW_PATH', 
                  'PROCESSED_PATH', 'ANALYSIS_PATH', 'RESULTS_PATH',
                  'ANALYSIS_DONE_PATH' )
