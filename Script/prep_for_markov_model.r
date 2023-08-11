# prep_for_markov_model.r

# function for setting different random seeds on each code section
set.seed.here.visible <- function( seed.base, seed.char )
{
  seed.add <- strtoi( substr( digest( seed.char, "xxhash32" ), 2, 8 ), 16 )
  seed.new <- seed.base + seed.add
  set.seed( seed.new )
  seed.new
}


# function for removing zeros from model data
model.x.remove.zeros <- function( x )
{
  x.pos.min <- min( x[ x > 0 ] )
  t( apply( x, 1, function( xr )
    if ( min( xr ) == 0 ) {
      xr <- xr + x.pos.min * 1e-3
      xr / sum( xr )
    }
    else
      xr
  ) )
}


# universal week alpha estimation
estimate.week.alpha.simple.n <- function( parabio.data, w, n.tis, pb )
{
  parabio.data.week.ids <- parabio.data[ parabio.data$week == w, "mouse" ]
  parabio.data.week.ids.in3 <- 
    names( table( parabio.data.week.ids ) == n.tis )[ 
      table( parabio.data.week.ids ) == n.tis ] %>% 
    as.numeric
  week.mouse <- parabio.data.week.ids.in3
  
  if ( w == 0 )                                                               
    week.popul <- gsub( pattern = "\\.(.*)\\.", 
                        replacement = ".celltype.", pb$pb.host.popul ) else                                        
      week.popul <- gsub( pattern = "\\.(.*)\\.", 
                          replacement = ".celltype.", pb$pb.hostdonor.popul )
  week.model.x <- lapply( pb$pb.tissue, function( tis ) {                   
    x <- t ( sapply( week.mouse, function( mou )                          
      as.vector( data.matrix( parabio.data[
        parabio.data$mouse == mou & parabio.data$tissue == tis,         
        week.popul                                                     
      ] ) )                                                           
    ) )
    rownames( x ) <- week.mouse
    colnames( x ) <- week.popul
    x <- sweep( x, 1, rowSums( x ), "/" )
    model.x.remove.zeros( x )                                              
  } )
  names( week.model.x ) <- pb$pb.tissue
  
  week.alpha <- t( sapply( pb$pb.tissue, function( tis )
    colMeans( week.model.x[[ tis ]] )
  ) )
  
  week.alpha
}
