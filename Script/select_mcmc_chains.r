# select_mcmc_chains.r

select_mcmc_chains <- function ( parabio.fit, mcmc.chain.n, type = "3mad" )
{
  fit.lp <- rep( NA, mcmc.chain.n )
  for ( i in 1 : mcmc.chain.n ) {
    try( fit.lp[ i ] <- get_posterior_mean( 
      parabio.fit )[ "lp__", i ] )
  }
  if ( type == "3mad" ) {
    fit.lp <- ifelse( fit.lp < 0, NA, fit.lp )
    if ( sum( !is.na( fit.lp ) ) == 0 ) {
      parabio.fit.chain <- c() } 
    if ( sum( !is.na( fit.lp ) ) == 1 ) {
      parabio.fit.chain <- which( !is.na( fit.lp ) )
    } else {
      parabio.fit.chain <- ( 1 : mcmc.chain.n )[
        which( abs( fit.lp - median( fit.lp, na.rm = TRUE) ) / 
                 stats::mad( fit.lp, na.rm = TRUE ) < 3 ) ]
    }
  } else if ( type == "max" ) {
    parabio.fit.chain <- which.max( fit.lp )  
  } else stop( "unknown type of sel.chains" )
  
  # print( fit.lp )
  # cat( "\nchain lp:", fit.lp, "\n"  )
  # cat( "chains selected:", parabio.fit.chain, "\n" )
  # cat( "chains rejected:", setdiff( 1 : mcmc.chain.n, parabio.fit.chain ),
  #      "\n\n" )
  
  return( parabio.fit.chain )
}
