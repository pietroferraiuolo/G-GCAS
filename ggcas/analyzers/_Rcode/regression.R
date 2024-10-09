library("minpack.lm")

regression <- function(data, method) {
  N <- seq(min(data), max(data), length.out=ceiling(1.5*sqrt(length(data))))
  hist_data <- hist(data, breaks = N, plot = FALSE)
  x <- hist_data$mids  # The midpoint of histogram bins
  y <- hist_data$counts  # The count in each bin
  if (method == "gaussian") {
    # Define Gaussian function
    gaussian <- function(x, a, mean, sd) {
      a * exp(-(x - mean)^2 / (2 * sd^2))
    }
    # Fit the Gaussian function to the histogram data
    fit_gaussian <- nlsLM(y ~ gaussian(x, a, mean, sd),
                          start = list(a = max(y), mean = mean(data), sd = sd(data)))
    out <- predict(fit_gaussian)
    coefficients <- coef(fit_gaussian)
  } else if (method == "boltzmann") {
    # Define the Boltzmann function
    boltzmann <- function(x, A1, A2, x0, dx) {
      (A1 - A2) / (1 + exp((x - x0) / dx)) + A2
    }
    fit_boltzmann <- nlsLM(y ~ boltzmann(x, A1, A2, x0, dx),
                           start = list(A1 = max(y), A2 = min(y), x0 = 0, dx = 1))
    out <- predict(fit_boltzmann)
    coefficients <- coef(fit_boltzmann)
  } else if (method == "exponential") {
    exponential <- function(x, a, b) {
      a * exp(-b * x)
    }
    fit_exponential <- nlsLM(y ~ exponential(x, a, b),
                             start = list(a = 10, b = 1))
    out <- predict(fit_exponential)
    coefficients <- coef(fit_exponential)
  } else {
    stop("Unknown method")
  }
  return(list(x=x, y=out, coeffs=coefficients))
}