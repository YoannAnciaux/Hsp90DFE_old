#' MLE on the DFE for FGM
#'
#' Maximum likelihood estimation (mle) for the parameters of the DFE density for
#' different form of the FGM.
#'
#' @inheritParams dfgm_martin
#' @param model String. "Martin" performs the mle on the density of the DFE from
#' Martin & Lenormand (2015) (see \code{\link{dfgm_martin}}). "Tenaillon" performs
#' the mle on the density of the DFE from Tenaillon (2014) (see \code{\link{dfgm_tenaillon}}).
#' Default "Martin".
#' @param start Numeric vector (named). Used as starting values in \code{\link{maxLik}{maxLik}}.
#' It must contains the parameters for the chosen \code{model}:
#' \itemize{
#'   \item{"Martin"}{ c(n = #, lambda = #, so = #)}
#'   \item{"Tenaillon"}{ c(n = #, lambda = #, so = #, alpha = #, Q = #)}
#' }
#' @param method String. Maximisation method. See \code{\link{maxLik}{maxLik}}.
#' Default "NM".
#' @param ... See \code{\link{maxLik}{maxLik}} all the possible parameters
#'
#' @return Returns the output of \code{\link[maxLik]{maxLik}} for the chosen model.
#'
#' @seealso \code{\link{dfgm_martin}}, \code{\link{dfgm_tenaillon}},
#' \code{\link{maxLik}{maxLik}}.
#'
#' @examples
#' #### model : "Martin" ####
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
#' res <- mle_dfgm(s, model = "Martin", start = initial_parameter,
#'                 method = "NM", constraints = list(ineqA = consA, ineqB = consB))
#'
#' #### model : "Tenaillon" ####
#' # Parameters
#' nsample <- 1000; n <- 4; lambda <- 0.5; so <- 1; alpha <- 1/2; Q <- 2;
#' # Simulated data
#' s <- rfgm(nsample, n, lambda, so, alpha, Q)
#' # Noise on initial parameters
#' initial_parameter <- c(n = ceiling(abs(rnorm(1, n))),
#'                        lambda = abs(rnorm(1, lambda, sd = 0.1)),
#'                        so = abs(rnorm(1, n)),
#'                        alpha = alpha,
#'                        Q = abs(rnorm(1, Q)))
#' # Constraints on parameters (not applioed for all methods see maxLik in maxLik package)
#' consA <- rbind(c(1, 0, 0, 0, 0),
#'                c(0, 1, 0, 0, 0),
#'                c(0, 0, 1, 0, 0),
#'                c(0, 0, 0, 1, 0),
#'                c(0, 0, 0, 0, 1))
#' consB <- c(0, 0, -max(s), 0, 0)
#' # MLE
#' res <- mle_dfgm(s, model = "Tenaillon", start = initial_parameter,
#'                 method = "NM", constraints = list(ineqA = consA, ineqB = consB))
#'
#' @source
#' \itemize{
#'   \item Tenaillon, O. (2014). The utility of Fisher's geometric model in
#'   evolutionary genetics. Annual review of ecology, evolution, and systematics,
#'   45, 179-201.
#'   \item Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#'   across environments: Fisher's geometrical model with multiple optima.
#'   Evolution, 69(6), 1433-1447.
#'  }
#'
#' @export
mle_dfgm <- function(s, model = "Martin", start, method = "NM", ...){

  coll <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(s, finite = T, any.missing = F, null.ok = F, add = coll)
  checkmate::assert_subset(model, choices = c("Martin", "Tenaillon"), empty.ok = F,
                           add = coll)
  checkmate::assert_vector(start, strict = T, any.missing = F, min.len = 3, max.len = 5,
               unique = F, names = "named", null.ok = F, add = coll)
  checkmate::reportAssertions(coll)

  # Log-Likelihood for the model "Tenaillon"
  if (model == "Martin") {
    checkmate::assert_names(names(start), permutation.of = c("n", "lambda", "so"))
    start <- start[c("n", "lambda", "so")]
    loglik <- function(param) {
      n <- param[1]
      lambda <- param[2]
      so <- param[3]
      ll <- sum(log(dfgm_martin(s, n, lambda, so)))
    }
  }
  # Log-Likelihood for the model "Tenaillon"
  if (model == "Tenaillon") {
    checkmate::assert_names(names(start), permutation.of = c("n", "lambda", "so", "alpha", "Q"))
    start <- start[c("n", "lambda", "so", "alpha", "Q")]
    loglik <- function(param) {
      n <- param[1]
      lambda <- param[2]
      so <- param[3]
      alpha <- param[4]
      Q <- param[5]
      ll <- sum(log(dfgm_tenaillon(s, n, lambda, so, alpha, Q)))
    }
  }

  # MLE
  maxLik::maxLik(loglik, start = start, method = method, ...)
}

#' Random generation of selection coefficient under FGM
#'
#' Generates selection coefficients of random mutations. First, the phenotype are
#' drawn from a multinormal distribution with means : 0, variances : lambda, and
#' covariances : 0. Second the coefficient of selections of these phenotypes are
#' computed using Fisher's Geometric Model (see Tenaillon (2014), Martin & Lenormand
#' (2015)).
#'
#' @param nsample Real number. Number of mutations to draw.
#' @inheritParams dfgm_tenaillon
#'
#' @return Vector of selection coefficients of length \code{nsample}.
#'
#' @examples
#' rfgm(10, 2, 0.1, 1, 1/2, 2)
#'
#' @source
#' \itemize{
#'   \item Tenaillon, O. (2014). The utility of Fisher's geometric model in
#'   evolutionary genetics. Annual review of ecology, evolution, and systematics,
#'   45, 179-201.
#'   \item Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#'   across environments: Fisher's geometrical model with multiple optima.
#'   Evolution, 69(6), 1433-1447.
#'  }
#'
#' @export
rfgm <- function(nsample, n, lambda, so, alpha, Q) {

  # Arguments
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_count(nsample, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_count(n, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_number(lambda, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(so, na.ok = F, finite = T, null.ok = F,
                           add = coll)
  checkmate::assert_number(alpha, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(Q, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::reportAssertions(coll)

  # Selection coefficients
  pheno_optium <- matrix(c((so / alpha)^(1/Q), numeric(n - 1)), nrow = nsample,
                         ncol = n, byrow = T)
  pheno_mutation <- MASS::mvrnorm(n = nsample, mu = numeric(n), Sigma = lambda * diag(1, n, n))
  so - alpha * apply(pheno_optium - pheno_mutation, 1, function(x) sqrt(sum(x^2))^Q)
}
