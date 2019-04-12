library(dplyr)

factor_cdf <- function(f, name = "Generic Factors", interval = 0.0001) {
    x <- seq(0.0, 1.0, interval)
    y <- mapply(f, x)
    factor_pdf <- cumsum(y / sum(y))
    return(data.frame(x = x, y = factor_pdf, name = name))
}

cubic_pdf <- function(x) {
    return((64 * (0.5 - x) ** 6) + (dnorm(x, 0.35, 0.04) / 10) +
           (dnorm(x, 0.65, 0.04) / 10))
}

quadratic_pdf <- function(x) {
    return((64 * (0.5 - x) ** 6) + (dnorm(x, 0.5, 0.04) / 10))
}

linear_pdf <- function(x) {
    return((64 * (0.5 - x) ** 6))
}

uniform_pdf <- function(x) {
    return(1.0)
}

write.csv(bind_rows(factor_cdf(linear_pdf, name = "linear"),
                    factor_cdf(quadratic_pdf, name = "quadratic"),
                    factor_cdf(cubic_pdf, name = "cubic"),
                    factor_cdf(uniform_pdf, name = "uniform")),
          "factor_cdfs.csv",
          row.names = FALSE)
