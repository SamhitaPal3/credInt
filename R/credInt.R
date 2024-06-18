#' Credible Interval for intended coverage
#'
#' @param X Design matrix
#' @param y Response variable
#' @param alpha Intended coverage
#' @param num_mcmc Number of posterior samples (default = 2000)
#'
#' @return Credible interval for all regression covariates and posterior samples
#' @export
#'
#' @examples credInt(X,y,alpha = 0.05,num_mcmc = 5000)
#'
#' @author Samhita Pal
credInt <- function(X,y,alpha,num_mcmc = 2000)
{
  n           <- nrow(X)
  p           <- ncol(X)
  l           <- glmnet::cv.glmnet(X,y,intercept = F)$lambda.min
  a_n         <- 2

  fit_new     <- matrix(0,ncol = num_mcmc, nrow = p)
  for(i in 1:num_mcmc)
  {
    sig2        <- 1/rgamma(1,2,2)
    temp        <- solve(t(X)%*%X + a_n*diag(p))
    post_sig    <- sig2*temp # Posterior Variance
    post_mean   <- temp%*%(t(X)%*%y)     # Posterior Mean

    post_samp   <- MASS::mvrnorm(1, mu = post_mean,
                                 Sigma = post_sig) # Full/dense
    new_resp    <- X%*%post_samp                # p x M

    fit_new[,i] <- as.vector(coef(glmnet::glmnet(X, new_resp, alpha = 1,
                                                     lambda = l,intercept = F,
                                                     type.measure = "mse"))[-1]) # Mapping
    print(i)
  }

  gamma       <- solve_gamma(l,alpha)

  ci.low      <- sapply(1:p, function(j) {quantile(fit_new[j,],prob = gamma/2)})
  ci.up       <- sapply(1:p, function(j) {quantile(fit_new[j,],prob = 1-gamma/2)})

  return(list("credInt" = cbind(ci.low,ci.up),"postDist" = fit_new))
}

