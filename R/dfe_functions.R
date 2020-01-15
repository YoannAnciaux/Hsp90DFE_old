#' DFE of FGM from Martin & Lenormand (2015).
#'
#' Density of the distribution of fitness effects of random mutations from \href{https://doi.org/10.1111/evo.12671}{Martin
#' & Lenormand (2015)}. The density corresponds to equation (2) in the article:
#' \deqn{f_s(s, n, \lambda, s_o) = \frac{2}{\lambda}%
#' f_{\chi_{n}^2} (\frac{2 (s_o-s)}{\lambda}, \frac{2 s_o}{\lambda})}{%
#' fs(s, n, \lambda, so) = 2/\lambda d\chi[2,n](2 (so - s) / \lambda, 2 so / \lambda)}
#' The density is encoded using the functions from the package \code{\link[distr]{distr}}
#'
#' @param s Vector of real numbers. Quantiles of the density. The density is 0
#' if \eqn{s \ge so}.
#' @param n Natural number. The dimensionality, i.e. the number of phenotypic
#' dimensions (traits) under selection.
#' @param lambda Positive real number. Mutational variance per trait among random
#' mutations scaled by the strength of selection. It scales the mean fitness effect
#' of mutations.
#' @param so Positive real number. Fitness distance between the wild-type (phenotype
#' without mutations) and the optimum.
#'
#' @return The density for each quantiles (selection coefficients) in a vector
#' format. The length of the vector is equal to the length of \code{s}.
#'
#' @examples
#' s <- seq(-2, 1, 0.1)
#' fs <- dfgm_martin(s, 2, 0.1, 1)
#'
#' @source Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#' across environments: Fisher's geometrical model with multiple optima.
#' Evolution, 69(6), 1433-1447.
#'
#'@export
dfgm_martin <- function(s, n, lambda, so) {

  #### Arguments checks ####
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(s, finite = T, any.missing = F, null.ok = F,
                            add = coll)
  checkmate::assert_count(n, na.ok = F, positive = T, null.ok = F, add = coll)
  checkmate::assert_number(lambda, na.ok = F, lower = .Machine$double.eps,
                           finite = T, null.ok = F, add = coll)
  checkmate::assert_number(so, na.ok = F, finite = T, null.ok = F,
                           add = coll)
  checkmate::reportAssertions(coll)

  #### Density ####
  C <- distr::Chisq(df = n, ncp = 2 * so / lambda)
  S <- so - lambda / 2 * C
  return(distr::d(S)(s))
}

#' DFE of FGM from Tenaillon (2014).
#'
#' Density of the distribution of fitness effects of random mutations from
#' \href{https://doi.org/10.1146/annurev-ecolsys-120213-091846}{Tenaillon (2014)}
#' adapted to show the same parameters as Martin & Lenormand (2015) plus the
#' parameters \code{Q} and \code{alpha}.
#' The density corresponds to the first equation in the paragraph *"Distribution of
#' mutation effects"* in the article.
#' The density is encoded using the functions from the package \code{\link[gsl]{gsl}}
#'
#' @inheritParams dfgm_martin
#' @param alpha Positive real number. Scaling parameter in the formulation of
#' Tenaillon et al. (2007). (Identifiability problem with \code{lambda})
#' @param Q Positive real number. Parameter that modifies the concavity of the
#' fitness decline from the optimum.
#'
#' @return The density for each quantiles (selection coefficients) in a vector
#' format. The length of the vector is equal to the length of \code{s}.
#'
#' @examples
#' s <- sample(seq(-2, 1.5, 0.1))
#' fs <- dfgm_tenaillon(s, 2, 0.1, 1, 1/2, 2)
#'
#' @source
#' \itemize{
#'   \item Tenaillon, O. (2014). The utility of Fisher's geometric model in
#'   evolutionary genetics. Annual review of ecology, evolution, and systematics,
#'   45, 179-201.
#'   \item Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#'   across environments: Fisher's geometrical model with multiple optima.
#'   Evolution, 69(6), 1433-1447.
#'   \item Tenaillon, O., Silander, O. K., Uzan, J. P., & Chao, L. (2007).
#'   Quantifying organismal complexity using a population genetic approach.
#'   PloS one, 2(2), e217.
#' }
#'
#'@export
dfgm_tenaillon <- function(s, n, lambda, so, alpha, Q) {

  #### Arguments checks ####
  coll <- checkmate::makeAssertCollection()
  checkmate::assert_numeric(s, finite = T, any.missing = F, null.ok = F,
                            add = coll)
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

  #### Density ####
  # The density is 0 is s > so, however the function gsl::bessel_In produces
  # NA/Nan/Inf when used with s > so. To overcome this numerical issue we compute
  # the density only for the values of s <= so and set the density to 0 for the
  # other.S
  d <- numeric(length(s))
  idx_s_infeqso <- which(s <= so)
  s <- s[idx_s_infeqso]

  # The formula from Tenaillon (2014) is separated in the for multiplicative terms
  # A B C D for readibility.
  A <- exp( -(((-s + so) / alpha)^(2/Q) + (so / alpha)^(2/Q)) / (2 * lambda) ) / (alpha * Q * lambda)
  B <- ((-s + so) / alpha)^((n/2 + 1 - Q)/Q)
  C <- (so / alpha)^((1 - n/2) / Q)
  D <- gsl::bessel_In(n = n/2 - 1, x = 1/lambda * ((-s + so) / alpha)^(1/Q) * (so / alpha)^(1/Q))
  d_s_infeqso <- A*B*C*D
  d[idx_s_infeqso] <- d_s_infeqso
  return(d)
}







