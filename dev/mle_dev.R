# File only for developpment of the functions
library(tidyverse)
nenv <- 8
lambda <- rexp(nenv, 5)
parameters_sim <- data.frame(nsample = rep(2000, nenv),
                             n = rgeom(nenv, 0.5) + 2,
                             lambda = lambda,
                             so = lambda * (rgeom(nenv, 0.2) + 2),
                             alpha = rep(1/2, nenv),
                             Q = rep(2, nenv))
parameters_sim <- asplit(parameters_sim, 1)
names(parameters_sim) <- as.character(1:nenv)

simulations <- as_tibble(lapply(parameters_sim,
                                function(e) as.vector(do.call("rfgm", as.list(e))))) %>%
  pivot_longer(everything(), names_to = "environment", values_to = "s")
data <- simulations
# e <- data %>% filter(environment == 1)
# s <- e$s
# # Starting parameters
# initial_parameter <- c(n = 2,
#                        lambda = 2 * abs(mean(s)) / 2,
#                        so = 2 * max(s))
# # Constraints on parameters
# consA <- rbind(c(1, 0, 0),
#                c(0, 1, 0),
#                c(0, 0, 1))
# consB <- c(0, 0, -max(s))
# # MLE
# res <- mle_dfgm(s, model = "Martin", start = initial_parameter,
#                 method = "NM", constraints = list(ineqA = consA, ineqB = consB))
# grid <- seq(min(s), max(s), length = 100)
# sample_dist_estim <- dplyr::tibble(s = grid,
#                                    density = dfgm_martin(grid, res$estimate[1], res$estimate[2], res$estimate[3]))
#
# ggplot() +
#   geom_density(aes(x = s), data = e) +
#   geom_line(aes(x = s, y = density), data = sample_dist_estim, color = "red") +
#   xlab("selection coefficient")

mle <- data %>%
  group_by(environment) %>%
  group_split() %>%
  mclapply(X = .,
           FUN = function(e) {
             s <- e$s
             # Starting parameters
             initial_parameter <- c(n = 2,
                                    lambda = 2 * abs(mean(s)) / 2,
                                    so = 2 * max(s))
             # Constraints on parameters
             consA <- rbind(c(1, 0, 0),
                            c(0, 1, 0),
                            c(0, 0, 1))
             consB <- c(0, 0, -max(s))
             # MLE
             res <- mle_dfgm(s, model = "Martin", start = initial_parameter,
                             method = "NM", constraints = list(ineqA = consA, ineqB = consB))
             grid <- seq(min(s), max(s), length = 100)
             list(mle = res,
                  dist_estim = tibble(environment = rep(e$environment[1], length(grid)),
                                      s = grid,
                                      density = dfgm_martin(grid, res$estimate[1], res$estimate[2], res$estimate[3])))
           }, mc.cores = 7)


tbl_param <- sapply(mle, FUN = function(e) {e$mle}, simplify = F)
tbl_dist_estim <- sapply(mle, FUN = function(e) {e$dist_estim}, simplify = F) %>%
  bind_rows()

ggplot(data, aes(x = s)) +
  geom_density() +
  geom_line(aes(x = s, y = density), data = tbl_dist_estim, color = "red") +
  facet_wrap(~ environment, scales = "free")

# # Sample the density with estimated parameters to be used in a plot
# grid <- seq(min(s), max(s), length = 100)
# if (model == "Martin") {
#   sample_dist_estim <- dplyr::tibble(s = grid,
#                                      density = dfgm_martin(grid, ml$estimate[1], ml$estimate[2], ml$estimate[3]))
# }
# if (model == "Tenaillon") {
#   sample_dist_estim <- dplyr::tibble(s = grid,
#                                      density = dfgm_tenaillon(grid, ml$estimate[1], ml$estimate[2], ml$estimate[3], ml$estimate[4], ml$estimate[5]))
# }
#
# # Output
# list(mle = ml,
#      dist_estim = sample_dist_estim)
#
# #TODO Add a class which can be used in a plot function
