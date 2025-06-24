#' Estimate Doublet Classification Using Central Matching and FDR Thresholding
#'
#' This function performs doublet detection on a dataset by applying Box-Cox transformation,
#' estimating null distributions, computing False Discovery Rate (FDR), and classifying
#' cells as "Doublet" or "Singlet" based on a specified threshold.
#'
#' @param dat A data frame containing at least the following columns:
#'   - `obs`: Observed values to be analyzed.
#'   - `label`: Ground truth labels ("Singlet" or "Doublet").
#'   - `barcode`: Unique cell identifiers.
#' @param truncation_point Numeric. The threshold below which cells are automatically labeled as "Singlet". Default is `0`.
#' @param pct0 Numeric vector of length 2. The lower and upper quantiles used for null distribution estimation. Default is `c(0.2, 0.6)`.
#' @param nulltype Integer. The type of null estimation model to use:
#'   - `2`: Symmetric model.
#'   - `3`: Asymmetric model allowing different variance estimates for left and right tails.
#' @param thres Numeric. FDR threshold for doublet classification. Default is `0.2`.
#'
#' @return A list containing:
#'   - `result`: A data frame with the original data and additional columns:
#'       - `box.cox.obs.truncated`: Transformed observations.
#'       - `cm.FDR.truncated`: Computed FDR values.
#'       - `cm.label.truncated`: Predicted classification ("Doublet" or "Singlet").
#'       - `cm.lfdr.truncated`: Local FDR estimates.
#'   - Additional parameters used in the estimation.
#'
#' @import dplyr
#' @import readr
#' @import stats
#' @import ggplot2
#' @import utils
#' @importFrom MASS boxcox
#' @importFrom splines ns
#'
#' @export

doublet_cm <- function(dat, truncation_point = 0, pct0 = c(0.2, 0.6), nulltype = 2, thres = 0.2) {

  # Ensure 'obs' exists in dataset
  if (!"obs" %in% colnames(dat)) stop("Column 'obs' not found in input data.")

  dat$filter <- dat$obs > truncation_point
  data.filtered <- dat %>% dplyr::filter(filter == TRUE)

  # Ensure filtered data contains 'obs'
  if (nrow(data.filtered) == 0) stop("Filtered data is empty after truncation. Adjust truncation_point.")

  # Apply Box-Cox transformation
  x <- data.filtered$obs
  res <- MASS::boxcox(lm(x ~ 1, y = TRUE))
  lambda <- res$x[which.max(res$y)]
  newx <- (x^lambda - 1) / lambda
  data.filtered$new.obs <- newx

  # Compute null estimation using an external function (assumed `cm_null_estimate`)
  est <- cm_null_estimate(zz = newx, bre = 120, df = 7, pct0 = pct0, nulltype = nulltype)

  # Compute Empirical CDF
  ecdf_function <- ecdf(newx)
  ecdf_values <- ecdf_function(newx)
  result <- data.frame(Observation = newx, ECDF = ecdf_values)

  # Compute FDR based on nulltype
  if (nulltype == 2) {
    mu0 <- est$delta.hat
    sigma <- est$sigma.hat
    pi0.hat <- est$p0

    ## ---------- goodness-of-fit on the CENTRAL window ----------
    lo.win <- quantile(newx, pct0[1])
    hi.win <- quantile(newx, pct0[2])
    window_vals <- newx[newx > lo.win & newx < hi.win]

    ## log-likelihood, AIC, BIC
    logLik_gaus <- sum(dnorm(window_vals,
                             mean = mu0, sd = sigma,
                             log  = TRUE))
    k   <- 2
    n   <- length(window_vals)
    AIC_val <- 2 * k - 2 * logLik_gaus
    BIC_val <- k * log(n) - 2 * logLik_gaus

    gof_list <- list(
      AIC     = AIC_val,
      BIC     = BIC_val
    )

    cat("AIC:", AIC_val, "\n")
    cat("BIC:", BIC_val, "\n")

    ## ------------------------------------------------------------

    p_x_greater <- 1 - pnorm(newx, mean = mu0, sd = sigma)
    # Visualization of histogram and null fit
    par(
      mar = c(5, 5, 4, 2),
      cex.lab = 1.4,
      cex.axis = 1.2,
      cex.main = 1.5,
      lwd = 2
    )

    hist(newx,
         breaks = 50,
         probability = TRUE,
         col = "gray85",
         border = NA,
         main = "Observed Distribution with Fitted Null (Singlet) Component",
         xlab = "Box-Cox Transformed Value",
         ylab = "Density",
         xlim = range(newx))

    x_vals <- seq(min(newx), max(newx), length.out = 1000)
    gauss_density <- dnorm(x_vals, mean = mu0, sd = ifelse(nulltype == 2, sigma, sigma.right))
    lines(x_vals, gauss_density, col = "darkblue", lwd = 3)
    grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
    legend("topright", legend = "Estimated Null (Singlet) Distribution",
           col = "darkblue", lwd = 3, bty = "n", cex = 1.1)
    abline(v = mu0, col = "#EE553D", lty = 2, lwd = 2)

  } else {
    mu0 <- est$delta.hat
    sigma.left <- est$sigma.left
    sigma.right <- est$sigma.right
    pi0.hat <- est$p0
    p_x_greater <- ifelse(newx > mu0,
                          1 - pnorm(newx, mean = mu0, sd = sigma.right),
                          1 - pnorm(newx, mean = mu0, sd = sigma.left))
    hist(newx,
         breaks = 50,
         probability = TRUE,
         col = "gray85",
         border = NA,
         main = "Observed Distribution with Asymmetric Null Fit",
         xlab = "Box-Cox Transformed Value",
         ylab = "Density",
         xlim = range(newx))

    # Sequence for left and right sides
    x_left <- seq(min(newx), mu0, length.out = 500)
    x_right <- seq(mu0, max(newx), length.out = 500)

    # Asymmetric Gaussian curves
    y_left <- dnorm(x_left, mean = mu0, sd = sigma.left)
    y_right <- dnorm(x_right, mean = mu0, sd = sigma.right)

    # Add the two density curves
    lines(x_left, y_left, col = "darkred", lwd = 3)
    lines(x_right, y_right, col = "darkred", lwd = 3)

    # Add grid and legend
    grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
    legend("topright", legend = "Estimated Asymmetric Null (Singlet)", col = "darkred", lwd = 3, bty = "n")

    # Vertical line at mu0
    abline(v = mu0, col = "blue", lty = 2, lwd = 2)
  }

  pi0.hat <- pmin(1,pi0.hat)

  # Compute FDR
  result$singlet.cdf.complement <- p_x_greater
  result$FDR <- (result$singlet.cdf.complement * pi0.hat) / (1 - result$ECDF)
  result$FDR[is.na(result$FDR) | is.infinite(result$FDR)] <- 0

  # Store results
  data.filtered$box.cox.obs.truncated <- result$Observation
  data.filtered$cm.FDR.truncated <- result$FDR
  data.filtered$cm.label.truncated <- ifelse(result$FDR < thres, "Doublet", "Singlet")

  # Merge results back to original dataset
  res.cm <- data.filtered %>%
    dplyr::select(barcode, obs, box.cox.obs.truncated, cm.FDR.truncated, cm.label.truncated)

  res <- merge(dat, res.cm, by = c("barcode", "obs"), all.x = TRUE)

  # Assign Singlet to truncated values
  res$cm.label.truncated[is.na(res$cm.label.truncated)] <- "Singlet"
  res$cm.FDR.truncated[is.na(res$cm.FDR.truncated)] <- max(res$cm.FDR.truncated, na.rm = TRUE)

  res <- res[, !(colnames(res) == "filter")]

  #lfdr recovery
  help <- res %>% arrange(cm.FDR.truncated)

  FDR <- help$cm.FDR.truncated
  cFDR <- seq_along(FDR) * FDR
  lfdr <- c(FDR[1], sapply(2:length(FDR), function(x) cFDR[x] - cFDR[x - 1]))

  lfdr_df <- data.frame(barcode = help$barcode, cm.lfdr.truncated = lfdr)
  res <- res %>% left_join(lfdr_df, by = "barcode")
  res$cm.lfdr.truncated <- ifelse(res$cm.lfdr.truncated>1, 1, res$cm.lfdr.truncated)


  # Corrected pi0 estimate
  pi0.correct <- (pi0.hat * length(newx) + sum(dat$obs <= truncation_point)) / length(dat$obs)

  # Count doublets and print the summary
  doublet_count <- sum(res$cm.label.truncated == "Doublet")
  total_cells <- nrow(res)
  doublet_proportion <- doublet_count / total_cells

  cat("Number of doublets called:", doublet_count, "\n")
  cat("Proportion of doublets:", round(doublet_proportion * 100, 2), "%\n")

  # Return final dataset
  if (nulltype == 3){

    return(list(result = res, pi0.hat.all = pi0.correct, pi0.hat = pi0.hat, mu0.hat = mu0, sigma.left = sigma.left,  sigma.right =  sigma.right))

  } else{

    return(list(result = res, gof = gof_list, pi0.hat.all = pi0.correct, pi0.hat = pi0.hat, mu0.hat = mu0, sigma.hat = sigma))

  }

}


#' Estimate the null distribution using central matching
#'
#' This function estimates the null distribution of a dataset by performing truncation,
#' Poisson regression, and either a symmetric or asymmetric quadratic approximation.
#'
#' @param zz Numeric vector. The dataset to be analyzed.
#' @param bre Integer. Number of breaks for histogram binning.
#' @param df Integer. Degrees of freedom for natural splines in Poisson regression.
#' @param pct0 Numeric vector of length 2. Defines the lower and upper quantile thresholds for selecting data used in null estimation.
#' @param nulltype Integer. The type of null estimation model:
#'   - `2`: Symmetric model using quadratic fitting.
#'   - `3`: Asymmetric model allowing different variance estimates for left and right tails.
#'
#' @return A list containing:
#'   - `delta.hat`: Estimated mode of the null distribution (for nulltype = 2).
#'   - `sigma.hat`: Estimated standard deviation of the null (for nulltype = 2).
#'   - `sigma.left`, `sigma.right`: Estimated left and right standard deviations (for nulltype = 3).
#'   - `p0`: Estimated null proportion for normalization.
#'
#' @export

cm_null_estimate <- function(zz, bre = 120, df = 7, pct0 = c(0.2, 0.6), nulltype) {

  # Set truncation limits
  lo <- min(zz)
  up <- max(zz)

  # Truncate values within the selected range
  zzz <- pmax(pmin(zz, up), lo)

  # Compute histogram
  breaks <- seq(lo, up, length = bre)
  zh <- hist(zzz, breaks = breaks, plot = FALSE)
  x <- (breaks[-1] + breaks[-length(breaks)]) / 2
  yall <- y <- zh$counts
  K <- length(y)
  N <- length(zz)

  # Poisson regression fit using natural splines
  X <- cbind(1, splines::ns(x, df = df))
  f <- glm(y ~ splines::ns(x, df = df), family = poisson)$fitted.values
  l <- log(f)

  # Compute central matching estimation
  imax <- which.max(l)
  xmax <- x[imax]

  # Define quantile thresholds
  pctlo <- pct0[1]
  pctup <- pct0[2]

  # Select data for null distribution estimation
  lo0 <- quantile(zz, pctlo)
  hi0 <- quantile(zz, pctup)
  nx <- length(x)
  i0 <- which(x > lo0 & x < hi0)
  x0 <- x[i0]
  y0 <- l[i0]

  ## If nulltype = 3 (Asymmetric model)
  if (nulltype == 3) {
    X00 <- cbind((x0 - xmax)^2, pmax(x0 - xmax, 0)^2)
    lr <- lm(y0 ~ X00)
    co <- lr$coef

    # Check if quadratic fit is valid
    cm_error <- is.na(co[3]) || co[2] >= 0 || (co[2] + co[3] >= 0)

    if (cm_error) {
      stop("CM estimation failed: Consider using nulltype = 2.")
    }

    X0 <- cbind(1, (x - xmax)^2, pmax(x - xmax, 0)^2)
    sigs <- 1 / sqrt(-2 * c(co[2], co[2] + co[3]))

    delta <- xmax
    sigma.left <- sigs[1]
    sigma.right <- sigs[2]

    l0 <- as.vector(X0 %*% co)
    f0 <- exp(l0)
    p0 <- sum(f0) / sum(f)
    f0 <- f0 / p0

    result <- list(delta.hat = delta, sigma.left = sigma.left, sigma.right = sigma.right, p0 = p0)

  } else if (nulltype == 2) {
    ## If nulltype = 2 (Symmetric model)
    X00 <- cbind(x0 - xmax, (x0 - xmax)^2)
    lr <- lm(y0 ~ X00)
    co <- lr$coef
    X0 <- cbind(1, x - xmax, (x - xmax)^2)

    # Compute peak and standard deviation
    xmaxx <- -co[2] / (2 * co[3]) + xmax
    sighat <- 1 / sqrt(-2 * co[3])

    # Compute null distribution function
    l0 <- as.vector(X0 %*% co)
    f0 <- exp(l0)
    p0 <- sum(f0) / sum(f)
    f0 <- f0 / p0

    result <- list(delta.hat = xmaxx, sigma.hat = sighat, p0 = p0)
  } else {
    stop("Invalid nulltype. Choose 2 or 3.")
  }


  return(result)
}






