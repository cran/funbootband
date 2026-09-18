## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  fig.width = 6, fig.height = 4,
  message = FALSE, warning = FALSE
)

# Keep the vignette fast on CRAN by using tiny B
is_cran <- !identical(tolower(Sys.getenv("NOT_CRAN")), "true")
B_demo  <- if (is_cran) 25L else 500L # bootstrap reps
k.coef_demo <- 4L

set.seed(1)

## ----iid-sim------------------------------------------------------------------
library(funbootband)

set.seed(1)
T <- 101L
n <- 30L
x <- seq(0, 1, length.out = T)
mu_true <- 0.7 * sin(2 * pi * x) - 0.2 * cos(4 * pi * x)

generate_iid_curve <- function() {
  mu_true +
    rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.30) * sin(2 * pi * x) +
    rnorm(1, sd = 0.20) * cos(2 * pi * x) +
    rnorm(1, sd = 0.15) * sin(4 * pi * x)
}

Y <- replicate(n, generate_iid_curve())

## ----iid-sim-plot-------------------------------------------------------------
# Fit prediction and confidence bands
fit_pred <- band(Y, type = "prediction", alpha = 0.10,
                 iid = TRUE, B = B_demo, k.coef = k.coef_demo)
fit_conf <- band(Y, type = "confidence", alpha = 0.10,
                 iid = TRUE, B = B_demo, k.coef = k.coef_demo)

## ----iid-plot, fig.cap="Calculated prediction (blue) and confidence (gray) bands."----
ylim  <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)

plot(x, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Normalized time", ylab = "Value",
     main = "Simultaneous bands (i.i.d.)")

matlines(x, Y, col = grDevices::adjustcolor("gray40", 0.25), lty = 1)
polygon(c(x, rev(x)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x, rev(x)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("darkorange", alpha.f = 0.30), border = NA)
lines(x, fit_pred$mean, lwd = 2)
lines(x, mu_true, col = "red", lwd = 2, lty = 2)

## ----clustered, eval=TRUE-----------------------------------------------------
library(funbootband)

set.seed(2)
K_subject <- 12L
m <- rep(c(2L, 3L, 4L), length.out = K_subject)
id <- rep(seq_len(K_subject), m)

subject_effect <- sapply(seq_len(K_subject), function(i) {
  rnorm(1, sd = 0.35) +
    rnorm(1, sd = 0.30) * sin(2 * pi * x) +
    rnorm(1, sd = 0.20) * cos(2 * pi * x)
})

within_subject_effect <- function() {
  rnorm(1, sd = 0.18) * sin(4 * pi * x) +
    rnorm(1, sd = 0.12) * cos(4 * pi * x)
}

Y <- sapply(seq_along(id), function(j) {
  mu_true + subject_effect[, id[j]] + within_subject_effect()
})

trial <- ave(id, id, FUN = seq_along)
colnames(Y) <- paste0("subject", id, "_trial", trial)


# Fit prediction and confidence bands
fit_pred <- band(Y, type = "prediction", alpha = 0.10, iid = FALSE,
                 id = id, B = B_demo, k.coef = k.coef_demo)
fit_conf <- band(Y, type = "confidence", alpha = 0.10, iid = FALSE,
                 id = id, B = B_demo, k.coef = k.coef_demo)

## ----clustered-meta-----------------------------------------------------------
fit_pred$meta[c("target", "weighting", "bootstrap_unit", "n_clusters")]

## ----clustered-plot-----------------------------------------------------------
ylim   <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)

plot(x, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Normalized time", ylab = "Value",
     main = "Simultaneous bands (clustered)")

matlines(x, Y, col = grDevices::adjustcolor("gray40", 0.20), lty = 1)
polygon(c(x, rev(x)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x, rev(x)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("darkorange", alpha.f = 0.30), border = NA)
lines(x, fit_pred$mean, lwd = 2)
lines(x, mu_true, col = "red", lwd = 2, lty = 2)

## ----clustered-future-coverage, eval = !is_cran-------------------------------
# generate_new_subject_curve <- function() {
#   new_subject_effect <-
#     rnorm(1, sd = 0.35) +
#     rnorm(1, sd = 0.30) * sin(2 * pi * x) +
#     rnorm(1, sd = 0.20) * cos(2 * pi * x)
#   new_curve_effect <-
#     rnorm(1, sd = 0.18) * sin(4 * pi * x) +
#     rnorm(1, sd = 0.12) * cos(4 * pi * x)
#   mu_true + new_subject_effect + new_curve_effect
# }
# 
# future_curves <- replicate(500L, generate_new_subject_curve())
# covered <- apply(future_curves, 2L, function(curve) {
#   all(curve >= fit_pred$lower & curve <= fit_pred$upper)
# })
# mean(covered)

## ----kcoef-mse-insample, eval = !is_cran--------------------------------------
# # MSE vs k.coef for the i.i.d. example (uses Y from above)
# 
# fourier_basis <- function(T, K) {
#   t <- 0:(T - 1L)
#   denom <- T - 1L
#   if (K == 0L) return(cbind(1))
#   cbind(
#     1,
#     sapply(1:K, function(k) cos(2*pi*k*t/denom)),
#     sapply(1:K, function(k) sin(2*pi*k*t/denom))
#   )
# }
# 
# reconstruct <- function(B, Y) {
#   coef <- qr.coef(qr(B), Y) # solve for all curves at once
#   B %*% coef
# }
# 
# mse <- function(A, B) mean((A - B)^2)
# 
# Ks <- c(10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L, 90L, 99L)
# T  <- nrow(Y)
# 
# tab <- do.call(rbind, lapply(Ks, function(K){
#   B <- fourier_basis(T, K)
#   Yhat <- reconstruct(B, Y)
#   data.frame(k.coef = K,
#              mse = mse(Y, Yhat),
#              pve = 100 * (1 - sum((Y - Yhat)^2) / sum((Y - mean(Y))^2)))
# }))
# row.names(tab) <- NULL
# print(tab)
# 
# # Plot
# op <- par(mar = c(4,4,2,1))
# plot(tab$k.coef, tab$mse, type = "b", xlab = "k.coef", ylab = "MSE",
#      main = "Fourier reconstruction error vs. k.coef")

## ----session-info-------------------------------------------------------------
sessionInfo()

