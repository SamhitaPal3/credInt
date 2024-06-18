#' Title
#'
#' @param lambda The cross-validated penalty parameter
#' @param alpha The desired level of confidence
#'
#' @return The choice of gamma for credible intervals
#' @export
#'
#' @examples solve_gamma(lambda = 0.02, alpha = 0.05)
#'
#' @author Samhita Pal
solve_gamma <- function(lambda, alpha)
{
  f <- function(gamma)
  {
    z_gamma_2 <- qnorm(1 - gamma / 2)
    diff      <- pnorm(lambda / 2 + z_gamma_2) - pnorm(lambda / 2 - z_gamma_2)

    diff - (1 - alpha)
  }

  # Use uniroot to find the root of the function f
  gamma_sol   <- uniroot(f, c(0, 1))$root

  return(gamma_sol)
}
