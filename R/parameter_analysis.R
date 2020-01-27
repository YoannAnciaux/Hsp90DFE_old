#' Confidence interval for MLE output
#'
#' Approximation of the confidence interval for object of class "maxLix" (See
#' \code{\link{maxLik}{maxLik}}) such as output from \code{\link{mle_dfgm}}.
#' The confidence interval are computed using the asymptotic normal approximation.
#' (Implemented following https://www.stat.umn.edu/geyer/5931/mle/mle.pdf)
#'
#' @param mle Object of class "maxLik". (See \code{\link{maxLik}{maxLik}}).
#' @param conf_level Number between 0 and 1. Default 0.95. Confidence level for the confidence
#' interval.
#' @param bonf Logical. Default True. If True applies a Bonferroni correction to
#' computes the critical value used for the confidence interval which depends on
#' the number of parameter. This provides simulatenous confidence interval. If False,
#' does not apply any correction, so the confidence intervals are not simultaneous
#' for motre than one parameter.
#'
#' @return Returns a matrix with the estimate, the standard error, the upper bound
#'  and the lower bound of the confidence interval.
#'
#' @seealso \code{\link{mle_dfgm}}
#'
#' @examples
#' # Parameters
#' nsample <- 1000; n <- 4; lambda <- 0.5; so <- 1; alpha <- 1/2; Q <- 2;
#' # Simulated data
#' s <- rfgm(nsample, n, lambda, so, alpha, Q)
#' # Noise on initial parameters
#' initial_parameter <- c(n = ceiling(abs(rnorm(1, n))),
#'                        lambda = abs(rnorm(1, lambda, sd = 0.1)),
#'                        so = abs(rnorm(1, n)))
#' # Constraints on parameters (not applioed for all methods see maxLik in maxLik package)
#' consA <- rbind(c(1, 0, 0),
#'                c(0, 1, 0),
#'                c(0, 0, 1))
#' consB <- c(0, 0, -max(s))
#' # MLE
#' mle <- mle_dfgm(s, model = "Martin", start = initial_parameter,
#'                 method = "NM", constraints = list(ineqA = consA, ineqB = consB))
#'
#' CI <- estim_CI(mle, bonf = F)
#' CI_simultaneous <- estim_CI(mle, bonf = T)
#'
#' @export
estim_CI <- function(mle, conf_level = 0.95, bonf = T){
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_class(mle, "maxLik", null.ok = FALSE, add = coll)
  checkmate::assert_number(conf_level, na.ok = FALSE, lower = 0, upper = 1,
                           finite = T, null.ok = FALSE, add = coll)
  checkmate::assert_flag(bonf, na.ok = FALSE, null.ok = FALSE, add = coll)
  checkmate::reportAssertions(coll)

  # Computes the critical value depending on the confidence level and the bonferroni correction
  if (bonf) {
    crit <- qnorm((1 + conf_level)/2)
  } else {
    crit <- qnorm(1 - (1 - conf_level)/2/length(mle$estimate))
  }

  # Computes asymptotic variance covariance matrix of the MLE from the Hessian matrix
  fisher_inf <- -mle$hessian
  asymptotic_var_mat_MLE <- solve(fisher_inf)

  # Computes lower and upper bounds of the confidence interval for each parameter.
  mat <- sapply(1:length(mle$estimate),
                function(p) mle$estimate[p] * c(1, 0, 1, 1) + c(0, 1, -1, 1) * crit * sqrt(asymptotic_var_mat_MLE[p, p]))
  colnames(mat) <- names(mle$estimate)

  as_tibble(mat) %>%
    mutate(stat = c("estimate", "SE", "lower_CI", "upper_CI"))
}
