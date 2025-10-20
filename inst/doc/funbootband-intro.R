## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  fig.width = 6, fig.height = 4,
  message = FALSE, warning = FALSE
)

# Keep the vignette fast on CRAN by using tiny B
is_cran <- !identical(tolower(Sys.getenv("NOT_CRAN")), "true")
B_demo  <- if (is_cran) 25L else 1000L # bootstrap reps
k.coef_demo <- 50L

set.seed(1)

## ----iid-sim------------------------------------------------------------------
library(funbootband)

set.seed(1)
T <- 200
n <- 10
x <- seq(0, 1, length.out = T)

# Simulate smooth Gaussian-process-like curves of equal length
mu  <- 10 * sin(2 * pi * x)
ell <- 0.12; sig <- 3
Kmat <- outer(x, x, function(s, t) sig^2 * exp(-(s - t)^2 / (2 * ell^2)))
ev <- eigen(Kmat + 1e-8 * diag(T), symmetric = TRUE)
Z  <- matrix(rnorm(T * n), T, n)
Y  <- mu + ev$vectors %*% (sqrt(pmax(ev$values, 0)) * Z)
Y  <- Y + matrix(rnorm(T * n, sd = 0.2), T, n) # observation noise

## ----iid-sim-plot-------------------------------------------------------------
# Fit prediction and confidence bands
fit_pred <- band(Y, type = "prediction", alpha = 0.11, iid = TRUE, B = B_demo, k.coef = k.coef_demo)
fit_conf <- band(Y, type = "confidence", alpha = 0.11, iid = TRUE, B = B_demo, k.coef = k.coef_demo)

## ----iid-plot, fig.cap="Calculated prediction (blue) and confidence (gray) bands."----
x_idx <- seq_len(fit_pred$meta$T)
ylim  <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)

plot(x_idx, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Index (Time)", ylab = "Amplitude",
     main = "Simultaneous bands (i.i.d.)")

matlines(x_idx, Y, col = "gray70", lty = 1, lwd = 1)
polygon(c(x_idx, rev(x_idx)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x_idx, rev(x_idx)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("gray40", alpha.f = 0.3), border = NA)
lines(x_idx, fit_pred$mean, col = "black", lwd = 1)

## ----clustered, eval=TRUE-----------------------------------------------------
library(funbootband)

set.seed(2)
T <- 200
m <- c(5, 5)
x <- seq(0, 1, length.out = T)

# Cluster-specific means
mu <- list(
  function(z) 8 * sin(2 * pi * z),
  function(z) 8 * cos(2 * pi * z)
)

# Generate curves with smooth within-cluster variation
Bm <- cbind(sin(2 * pi * x), cos(2 * pi * x))
gen_curve <- function(k) {
  sc <- rnorm(ncol(Bm), sd = c(2.0, 1.5))
  mu[[k]](x) + as.vector(Bm %*% sc)
}

Ylist <- lapply(seq_along(m), function(k) {
  sapply(seq_len(m[k]), function(i) gen_curve(k) + rnorm(T, sd = 0.6))
})
Y <- do.call(cbind, Ylist)
colnames(Y) <- unlist(mapply(
  function(k, mk) paste0("C", k, "_", seq_len(mk)),
  seq_along(m), m
))


# Fit prediction and confidence bands
fit_pred <- band(Y, type = "prediction", alpha = 0.11, iid = FALSE, B = B_demo, k.coef = k.coef_demo)
fit_conf <- band(Y, type = "confidence", alpha = 0.11, iid = FALSE, B = B_demo, k.coef = k.coef_demo)

## ----clustered-plot-----------------------------------------------------------
x_idx <- seq_len(fit_pred$meta$T)
ylim   <- range(c(Y, fit_pred$lower, fit_pred$upper), finite = TRUE)

plot(x_idx, fit_pred$mean, type = "n", ylim = ylim,
     xlab = "Index (Time)", ylab = "Amplitude",
     main = "Simultaneous bands (clustered)")

matlines(x_idx, Y, col = "gray70", lty = 1, lwd = 1)
polygon(c(x_idx, rev(x_idx)), c(fit_pred$lower, rev(fit_pred$upper)),
        col = grDevices::adjustcolor("steelblue", alpha.f = 0.25), border = NA)
polygon(c(x_idx, rev(x_idx)), c(fit_conf$lower, rev(fit_conf$upper)),
        col = grDevices::adjustcolor("gray40", alpha.f = 0.3), border = NA)
lines(x_idx, fit_pred$mean, col = "black", lwd = 1)

## ----kcoef-mse-insample, eval = !is_cran--------------------------------------
# # MSE vs k.coef for the i.i.d. example (uses Y from above)
# 
# fourier_basis <- function(T, K) {
#   t <- seq_len(T)
#   if (K == 0L) return(cbind(1))
#   cbind(
#     1,
#     sapply(1:K, function(k) cos(2*pi*k*t/T)),
#     sapply(1:K, function(k) sin(2*pi*k*t/T))
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

