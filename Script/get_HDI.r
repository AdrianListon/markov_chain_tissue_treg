#  get_HDI.r

get_HDI_4_values <- function( MCMC.sample, credMass.big = 0.80, credMass.small = 0.5 ) {
  credMass <- NULL
  sorted.sample <- sort( MCMC.sample )
  ciIdxInc.b <- ceiling( credMass.big * length( sorted.sample ) )
  ciIdxInc.s <- ceiling( credMass.small * length( sorted.sample ) )
  nCIs.b <- length( sorted.sample ) - ciIdxInc.b
  nCIs.s <- length( sorted.sample ) - ciIdxInc.s
  ciWidth.b <- rep( 0, nCIs.b )
  ciWidth.s <- rep( 0, nCIs.s )
  for ( i in 1:nCIs.b ) {
    ciWidth.b[ i ] <- sorted.sample[ i + ciIdxInc.b ] - sorted.sample[ i ]
  }
  for ( i in 1:nCIs.s ) {
    ciWidth.s[ i ] <- sorted.sample[ i + ciIdxInc.s ] - sorted.sample[ i ]
  }
  HDImin.b <- sorted.sample[ which.min( ciWidth.b ) ]
  HDImax.b <- sorted.sample[ which.min( ciWidth.b ) + ciIdxInc.b ]
  HDImin.s <- sorted.sample[ which.min( ciWidth.s ) ]
  HDImax.s <- sorted.sample[ which.min( ciWidth.s ) + ciIdxInc.s ]
  HDIlim.4 <- sort( c( HDImin.b, HDImin.s, HDImax.s, HDImax.b ) )
  return( HDIlim.4 )
}
