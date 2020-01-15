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

