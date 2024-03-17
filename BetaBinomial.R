#install.packages("rmutil")
#install.packages("boot")

# calculate likelihoods for mixture EM

library("rmutil")
library("boot")

SumLog <- function(a, b)
{
    if (a > b)
         return( a + log( 1.0 + exp(b - a) ) )
    else
         return( b + log( 1.0 + exp(a - b) ) )
}

#MixBetaBinomialLogDensity - returns log density given parameters and observations

#R-function: dbetabinom(observed, N, m[1], v[1], log = True) - returns betabinomial density.

#observed = vector of observations (number of loF carriers per gene)
#N = number of sequenced samples
#m[1] = average probability of beta distribution
#v[1] = dispersion of beta distribution
#*pi = mixing proportions for the beta-binomial components


MixBetaBinomialLogDensity <- function( observed, model, N )
{
  #model parameters
    pi <- model[1,]
    alpha <- model[2,]
    beta <- model[3,]
  
  #mean,variance
    m <- alpha / (alpha + beta)
    v <- beta / (1 - m)

    lk <- log(pi[1]) + dbetabinom(observed, N, m[1], v[1], log = TRUE)

    if (length(pi) == 1) 
       return(lk);

    for (i in 2:length(pi))
       lk <- SumLog(lk, log(pi[i]) + dbetabinom(observed, N, m[i], v[i], log = TRUE))

    return (lk);
}


MixBetaBinomialLogLikelihood <- function( observed, model, N )
{
    u <- unique(observed)
    d <- sapply(u, MixBetaBinomialLogDensity, model, N)
    densities <- d[match(observed,u)]

    return (sum(densities))
}


MixBetaBinomialUpdate <- function( observed, model, N)
{
    pi <- model[1,]
    alpha <- model[2,]
    beta <- model[3,]

    m <- alpha / (alpha + beta)
    v <- beta / (1 - m)

    u <- unique(observed)
    l <- matrix(data = 0, nrow=length(observed), ncol=length(pi))

    for (i in 1:length(u))
       {
       sumllk <- 0

       for (j in 1:length(pi))
           {
           l[i,j] <- log(pi[j]) + dbetabinom(u[i], N, m[j], v[j], log = TRUE)
           sumllk <- ifelse(j == 1, l[i,1], SumLog(sumllk, l[i,j]))
           }
      
       l[i,] <- exp(l[i,] - sumllk)
       }

    newProportions <- rep(0, length(pi));
    moment_x <- rep(0, length(alpha));
    moment_xx <- rep(0, length(beta));

    for (i in 1:length(observed))
       {
       posterior <- l[match(observed[i],u),]

       newProportions <- newProportions + posterior
       moment_x <- moment_x + posterior * observed[i]
       moment_xx <- moment_xx + posterior * observed[i] * observed[i]
       }

    newProportions <- sapply(newProportions, max, 1e-30)
    moment_x <- sapply(moment_x, max, 1e-30)
    moment_xx <- sapply(moment_xx, max, 1e-30)

    pi <- newProportions / length(observed)

    moment_x <- moment_x / newProportions
    moment_xx <- moment_xx / newProportions

    q <- N * ( moment_xx / moment_x - moment_x - 1) + moment_x

    alpha <- (N * moment_x - moment_xx) / q     
    beta <- (N - moment_x) * (N - moment_xx / moment_x) / q
    

    return(rbind(pi, alpha, beta))
}

FindBetaBinomialModel <- function( observed, N, components, refineWithOptim = F, quiet = F, return_refined = F )
{
    proportions <- rep(1/components, components)
    alpha <- 1:components
    beta <- 1 + 1 / ( mean(observed) / N )

    model <- rbind(proportions, alpha, beta)

    llk0 <- MixBetaBinomialLogLikelihood(observed, model, N)
    iter <- 1

    while (1)
       {
       model <- MixBetaBinomialUpdate(observed, model, N)

       llk1 <- MixBetaBinomialLogLikelihood(observed, model, N)

       if ( abs(llk0 - llk1) < abs( (llk0 + llk1) * 1e-6) )
          {
          if (quiet == F) cat("Best Log Likelihood: ", llk1, "after ", iter, "iterations\n")
          if (refineWithOptim)
             {
             v <- optim(ModelToVector(model), OptimWrapper, method = "L-BFGS-B", observed = observed, N = N)
             model <- VectorToModel(v$par)
             if (quiet == F) cat("  with refinement: ", -v$value, "after ", v$counts[1], "additional iterations\n")
            }
          if(refineWithOptim & return_refined)
            {
            return(-v$value)
          }
#            return(llk1)
          return(model)
          }
       
       llk0 <- llk1
       iter <- iter + 1
       }
}

PointPrediction <- function( observed, sampleSize, model )
{ 
    pi <- model[1,]
    alpha <- model[2,]
    beta <- model[3,]

    m <- alpha / (alpha + beta)
    v <- beta / (1 - m)
    prediction <- 0

    for (i in 1:length(pi))
       prediction <- prediction + pi[i] * (1 - pbetabinom(observed - 0.5, sampleSize, m[i], v[i]))

    return (prediction);
}

MatrixPrediction <- function (observed, sampleSize, model )
{
    prediction <- matrix( nrow = length(observed), ncol = length(sampleSize), dimnames = list(observed, sampleSize))

    for (i in 1:length(observed) )
      for (j in 1:length(sampleSize) )
          prediction[i,j] <- PointPrediction( observed[i], sampleSize[j], model )

    return(prediction)
}

QuickPredictions <- function(model)
{
    MatrixPrediction( c(1, 10, 100, 1000), c(1500, 15000, 150000), model)
}

LoopThroughModels <- function ( observedData, sampleSize, maxSize = 10 )
{

    llk <- rep(0, maxSize)
 
    for (modelSize in 1:maxSize)
       {
       model <- FindBetaBinomialModel( observedData, sampleSize, modelSize, refineWithOptim = T, quiet = T)
       llk[modelSize] <- MixBetaBinomialLogLikelihood( observedData, model, sampleSize)

       parameters <- 3 * modelSize - 1
       bicPenalty <- log(length(observedData)) * parameters

       cat(modelSize, "components\t llk =", llk[modelSize], "\taic = ", -2 * llk[modelSize] + 2 * parameters, "\tbic = ", -2 * llk[modelSize] + bicPenalty, "\n")
       }

    return(llk)
}

PlotPredictions <- function(model, title = "LoF Accumulation Plot")
{
   x <- as.numeric(labels(model)[[2]])
   z <- as.numeric(labels(model)[[1]])

   nlines <- length(z)

   colors <- rainbow(nlines)

   plot(model[1,] ~ x, type = "b", xlab = "Sample Size", 
        ylab = "Genes with LoF", lwd = 3, ylim = c(0,1),
        col = colors[1], main = title)

   for (i in 2:nlines)
      lines(model[i,] ~ x, type = "b", lwd = 3, col = colors[i])

   grid();

   legend("bottomright", title = "LoF carriers", lwd = 3, col = colors, 
           legend = z)
}

ModelToVector <- function(model)
   {
   v <- rep(0, length(model[1,]) - 1)

   sum = 1.0
   if (length(v) > 0)
      for (i in 1:length(v))
          {
          v[i] <- logit(model[1,i] / sum)
          sum <- sum - model[1,i]
          }   

   v <- c(v, log(model[2,]))
   v <- c(v, log(model[3,]))

   return(v)
   }

VectorToModel <- function(vector)
   {
   dim <- (length(vector) + 1) / 3

   pi <- rep(0, dim)
   
   sum <- 1.0
 
   if (dim > 1)
      for (i in 1:(dim-1))
         {
         pi[i] <- inv.logit(vector[i]) * sum
         sum <- sum - pi[i]
         }
   pi[dim] = 1.0 - sum(pi)

   alpha <- exp(vector[(dim):(2*dim-1)])
   beta <- exp(vector[(2*dim):(3*dim-1)])

   return (rbind(pi,alpha,beta))
   }

OptimWrapper <- function(par, observed, N)
   {
   model <- VectorToModel(par)

   -MixBetaBinomialLogLikelihood(observed, model, N)
   }

PlotDataAndFit <- function(observed, model, N, zoom = 100)
   {
   breaks <- 0:(zoom+1)

   hist( sapply(observed, min, zoom), breaks = breaks - 0.5, col = "pink", 
         main = "Distribution of LoF", xlab = "LoF Count")

   breaks <- 0:zoom
   fitted <- exp(sapply(breaks, MixBetaBinomialLogDensity, model = model, N = N))
   fitted[length(breaks)] <- PointPrediction(zoom, N, model = model)

   lines( breaks, fitted * length(observed), lwd = 3, col = "blue")
   }

