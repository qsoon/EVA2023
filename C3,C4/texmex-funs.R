### functions modified for keef-et-al.R ###
 
### Unif margin -> Laplace margin ###
laplace.transform <- function(p) ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p)))
laplace.inv.cdf <- function(q) ifelse(q < 0, exp(q)/2, 1 - 0.5 * exp(-q))
gumbel.transform <- function(p) -log(-log(p))

lap2gumb <- function(data){
  res <- apply(data$data$transformed,2,function(xx)-log(-log(laplace.inv.cdf(xx))))
  colnames(res) <- colnames(data$data$simulated)
  res
}


### Codes below are from the `texmex` package ### 

mexTransform <-
  function(x, margins, r=NULL, method = "mixture", divisor = "n+1", na.rm = TRUE){ # if required, r output from makeReferenceMarginalDistribution
    
    # if (!is.element(method, c("mixture", "empirical")))
    #   stop("method should be either 'mixture' or 'empirical'")
    # if (!is.element(divisor, c("n", "n+1")))
    #   stop("divisor can be 'n' or 'n+1'")
    
    if (is.null(r)){
      r <- x
      r$transData <- lapply(1:dim(x$data)[2],function(i)x$data[,i])
    }
    
    transFun <- function(i, x, r, mod, th, divisor, method){
      x <- x[, i]
      r <- r[[i]]
      mod <- mod[[i]]
      th <- th[i]
      
      if (divisor == "n") divisor <- length(r)
      else if (divisor == "n+1") divisor <- length(r) + 1
      
      ox <- order(x)
      
      r <- sort(r)
      run <- rle(r)
      p <- cumsum(run$lengths) / divisor
      p <- rep(p, run$lengths) # this calculated from r
      
      Femp <- p[sapply(x,function(y) which.min(abs(r-y)))]
      if (method == "mixture"){
        sigma <- exp(mod$coefficients[1])
        xi <- mod$coefficients[2]
        
        Para <- (1 + xi * (x - th) / sigma) ^ (-1 / xi) # this calculated from model fitted to r but data values are x
        Para <- 1 - mean(r > th) * Para
        res <- ifelse(x <= th, Femp, Para)
      }
      else res <- Femp
      
      res[ox] <- sort(res)
      res
    } # Close transfun
    
    res <- sapply(1:ncol(x$data), transFun,
                  x = x$data, r = r$transData, mod = r$models, th = r$mth,
                  divisor = divisor, method=method)
    
    dimnames(res) <- list(NULL, names(r$models))
    
    # x$transformed <- margins$p2q(res)
    # modify the transformation stage, as we know the exact marginal distribution
    # Gumbel margin to Laplace margin:
    x$transformed <- apply(x$data,2,function(cols) laplace.transform(exp(-exp(-cols))))
    
    invisible(x)
  }

### replace the invisible function `texmex::mexTransform` by the function defined above ###
assignInNamespace("mexTransform", mexTransform, "texmex")


`predict.mex` <-
  function(object, which, pqu = .99, nsim = 1000, trace=10, smoothZdistribution=FALSE, ...){
    theCall <- match.call()
    
    # Class can be either mex or bootmex
    theClass <- class(object)
    theClass <- theClass[1] # Avoid R CMD check that doesn't like indexing class(x) at 1
    if (! theClass %in% c("mex", "bootmex")){
      stop("object must have class 'mex' or 'bootmex'")
    }
    
    if (theClass == "bootmex" ){
      which <- object$which
      migpd <- object$simpleMar
      margins <- object$margins
      constrain <- object$constrain
      dall <- mexDependence( migpd , which=which , dqu=object$dqu, margins = margins[[1]], constrain=constrain )
    } else {
      which <- object$dependence$which
      if(is.null(object$margins$referenceMargin)){
        migpd <- object$margins
      } else {
        migpd <- object$margins$referenceMargin
      }
      margins <- object$dependence$margins
      constrain <- object$dependence$constrain
      dall <- object
    }
    
    ################################################################
    MakeThrowData <- function(dco,z,coxi,coxmi,data){
      ui <- runif(nsim , min = max(c(migpd$mqu[which], pqu)))
      y <- margins$p2q(ui)
      distFun <- margins$q2p
      
      z <- as.matrix(z[ sample( 1:( dim( z )[ 1 ] ), size=nsim, replace=TRUE ) ,])
      if(smoothZdistribution){
        z <- apply(z,2,function(x)x + rnorm(length(x),0,bw.nrd(x)))
      }
      ymi <- sapply( 1:( dim( z )[[ 2 ]] ) , makeYsubMinusI, z=z, v=dco , y=y )
      
      # xmi <- apply( ymi, 2, distFun )
      # 
      # xi <- u2gpd( ui, p = 1 - migpd$mqu[which], th = migpd$mth[ which ],
      #              sigma = coxi[ 1 ], xi = coxi[ 2 ] )
      # 
      # for( i in 1:( dim( xmi )[[ 2 ]] ) ){
      #   xmi[, i ] <- revTransform( xmi[ ,i ], as.matrix(data[,-which])[, i ],
      #                              th = migpd$mth[ -which ][ i ],
      #                              qu = migpd$mqu[ -which ][ i ],
      #                              sigma=coxmi[ 1,i ], xi=coxmi[ 2,i ] )
      # }
      #### lines commented above are unnecessary since we know the exact marginal distribution
      xmi <- apply(ymi,2,function(xx)-log(-log(laplace.inv.cdf(xx))))
      xi <- -log(-log(laplace.inv.cdf(y)))
      
      sim <- data.frame( xi , xmi , y, ymi)
      names( sim ) <- c( colnames( migpd$data )[ which ], colnames( migpd$data )[ -which ],
                         paste0(c(colnames( migpd$data )[ which ], colnames( migpd$data )[ -which ]),".trans"))
      sim[,dim(sim)[2]+1] <- y > apply(ymi,1,max) # condlargest extra column
      sim
    }
    
    ################################################################
    makeYsubMinusI <- function( i, z, v , y ){
      v <- v[ , i ]
      z <- z[ , i ]
      if ( !is.na( v[ 1 ] ) ){
        if( v[ 1 ] < 10^(-5) & v[ 2 ] < 0 ){
          if( v[ 4 ] < 10^(-5 ) ) d <- 0
          else d <- v[ 4 ]
          a <- v[ 3 ] - d * log( y )
        }
        else a <- v[ 1 ] * y
      } # close if( !is.na...
      else a <- NA
      a + ( y^v[ 2 ] ) * z
    }
    
    ###############################################################
    if (theClass == "bootmex"){
      # The function lfun does most of the work
      lfun <- function( i , bo, pqu, nsim , migpd, which ){
        if ( i %% trace == 0 ) cat( i, "sets done\n" )
        
        res <- MakeThrowData(dco=bo[[ i ]]$dependence,z=bo[[ i ]]$Z, coxi = bo[[i]]$GPD[,which],
                             coxmi = as.matrix(bo[[ i ]]$GPD[,-which]),
                             data = bo[[i]]$Y)
        res <- res[,1:((dim(res)[2]-1)/2)]
        res
      }
      
      bootRes <- lapply( 1:length( object$boot ) , lfun ,
                         migpd=migpd, pqu=pqu, bo = object$boot, nsim=nsim,
                         which = which )
      # bootRes contains the bootstrap simulated complete vectors X on the original
      # scale of the data, conditional on having the _which_ component above the pqu quantile.
    } else {
      bootRes <- NULL
    }
    
    ##########################################################################
    # Get a sample using the point estimates of the parameters
    # that are suggested by the data
    
    cox <- coef(migpd)[3:4, which]
    coxmi <- as.matrix(coef(migpd)[3:4, -which])
    
    sim <- MakeThrowData(dco=dall$dependence$coefficients, z=dall$dependence$Z,
                         coxi=cox, coxmi=coxmi, data=migpd$data)
    CondLargest <- sim[,dim(sim)[2]]
    transformed <- sim[,(((dim(sim)[2]-1)/2)+1):(dim(sim)[2]-1)]
    sim <- sim[,1:((dim(sim)[2]-1)/2)]
    
    m <- 1 / ( 1 - pqu ) # Need to estimate pqu quantile
    zeta <- 1 - migpd$mqu[ which ] # Coles, page 81
    pth <- migpd$mth[ which ] + cox[ 1 ] / cox[ 2 ] * ( ( m*zeta )^cox[ 2 ] - 1 )
    
    data <- list( real = data.frame( migpd$data[, which ], migpd$data[, -which] ), simulated = sim, pth=pth,CondLargest=CondLargest, transformed = transformed)
    names(data$real)[1] <- colnames(migpd$data)[which]
    
    res <- list( call = theCall , replicates = bootRes, data = data,
                 which = which, pqu = pqu,
                 mth=c( migpd$mth[ which ], migpd$mth[ -which ] ),
                 gpd.coef = coef(migpd)[,c(which,(1:dim(data$real)[2])[-which])])
    
    oldClass( res ) <- "predict.mex"
    
    res
  }

u2gpd <-
  function(u, p = 1, th=0, sigma, xi) {
    qgpd((1 - u) / p, sigma=sigma, u=th, xi=xi, lower.tail=FALSE)
  }

revTransform <-
  function (x, data, qu, th = 0, sigma = 1, xi = 0,
            method = c("mixture", "empirical")) {
    method <- match.arg(method)
    
    n <- length(data)
    probs <- (1:n)/(n + 1)
    
    px <- vapply(x,
                 function(x, p) {
                   p[[which.min(abs(x-p))]]
                 }, 0, p=probs)
    
    px <- as.integer(round(px * (1 + n)))
    res <- sort(data)[px]
    
    if (method == "mixture"){
      # Real data contain ties which can cause x[res > th] < qu, res[res < th] > qu
      i.x <- x >= qu
      i.r <- res > th
      i.rx <- apply(cbind(i.x, i.r), 1, all)
      if (sum(i.rx > 0)){
        wh <- u2gpd(x[i.rx], p=1-qu, th=th, sigma=sigma, xi=xi)
        rth <- res[i.rx]
        o <- order(rth)
        rth <- rth[o]
        rth[length(rth):(length(rth) - length(wh) + 1)] <- rev(sort(wh))
        rth <- rth[order(o)]
        res[i.rx] <- rth
      }
    }
    
    res[order(x)] <- sort(res)
    res
  }
