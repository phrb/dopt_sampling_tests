library(Rcpp)
library(rsm)
library(dplyr)
library(AlgDesign)
library(uuid)
library(MASS)
library(logging)

sourceCpp("rcpp_biased_sample.cpp")
sourceCpp("rcpp_squared_exponential_kernel.cpp")

generate_design_codings <- function(design){
    design_codings <- list()

    for (i in 1:ncol(design)) {
        design_codings[[paste("x", i, sep = "")]] <- formula(paste(paste("x",
                                                                         i,
                                                                         "~",
                                                                         sep = ""),
                                                                   paste("(X",
                                                                         i,
                                                                         "-0.5)/0.5",
                                                                         sep = ""),
                                                                   sep = ""))
    }

    return(design_codings)
}

encode_design <- function(design){
    design_codings <- generate_design_codings(design)
    return(coded.data(design, formulas = design_codings))
}

sample_gaussian_process <- function(data, amplitude, lengthscale, ...) {
    experiments <- data %>%
        dplyr::select(...) %>%
        as.matrix()

    loginfo("> Computing covariance matrix")

    covariance_matrix <- rcpp_squared_exponential_kernel(experiments,
                                                         experiments,
                                                         amplitude,
                                                         lengthscale)

    loginfo("> Sampling from multivariate normal")

    data$Y <- unname(mvrnorm(n = 1,
                             rep(0.0, nrow(covariance_matrix)),
                             covariance_matrix))

    data
}

biased_sample <- function(sampling_cdf, size, factors, name){
    rcpp_biased_sample(sampling_cdf,
                       size,
                       factors) %>%
        data.frame() %>%
        mutate(name = name) %>%
        encode_design()
}

compare_sampling_strategies <- function (uniform_strategy_cdf,
                                         biased_strategy_cdf,
                                         factors,
                                         testing_sample_size,
                                         selection_sample_size,
                                         design_size,
                                         regression_formula,
                                         filter_thresholds,
                                         gaussian_amplitude,
                                         gaussian_lengthscale){
    loginfo("> Starting to generate samples")

    testing_sample <- biased_sample(uniform_strategy_cdf,
                                    testing_sample_size,
                                    factors,
                                    "testing")

    selection_sample <- biased_sample(biased_strategy_cdf,
                                      selection_sample_size,
                                      factors,
                                      "selection")

    loginfo("> Generating random design")

    uniform_design <- biased_sample(uniform_strategy_cdf,
                                    design_size,
                                    factors,
                                    "uniform")

    loginfo("> Generating biased design")

    biased_design <- optFederov(regression_formula,
                                selection_sample,
                                nTrials = design_size)$design %>%
                                                     mutate(name = "biased")

    loginfo("> Sampling from gaussian process")

    sampled_data <- sample_gaussian_process(rbind(testing_sample,
                                                  uniform_design,
                                                  biased_design),
                                            gaussian_amplitude,
                                            gaussian_lengthscale,
                                            -name)

    uniform_significant_factors <- run_anova_tests(sampled_data %>%
                                                   subset(name == "uniform"),
                                                   regression_formula)

    biased_significant_factors <- run_anova_tests(sampled_data %>%
                                                  subset(name == "biased"),
                                                  regression_formula)
}

run_experiments <- function(){
    logging::basicConfig(level = "DEBUG")

    cdf_data <- read.csv("factor_cdfs.csv",
                         header = TRUE)

    factors <- 30
    amplitude <- 2.0

    significance_probability <- 0.15
    significance_variability <- list(max = 30, min = 10)
    significance_noise <- 0.02
    lengthscale <- replicate(factors,
                             ifelse(runif(1) < (1.0 - significance_probability),
                                    runif(1,
                                          min = significance_variability$min,
                                          max = significance_variability$max),
                                    rnorm(1, mean = 0.2, sd = significance_noise)))

    uniform_strategy_cdf <- cdf_data %>%
        subset(name == "uniform") %>%
        dplyr::select(x, y) %>%
        as.matrix()

    biased_strategy_cdf <- cdf_data %>%
        subset(name == "linear") %>%
        dplyr::select(x, y) %>%
        as.matrix()

    factor_names <- paste("x", 1:factors, sep = "")
    model_formulas <- list(linear = factor_names %>%
                               paste(sep = "", collapse = " + ") %>%
                               paste("~ ", ., sep = "") %>%
                               formula(),
                           quadratic = factor_names %>%
                               paste(sep = "", collapse = " + ") %>%
                               paste("+ ") %>%
                               paste(., paste("I(",
                                              factor_names,
                                              "^2)",
                                              collapse = " + ",
                                              sep = "")) %>%
                               paste("~ ", .) %>%
                               formula())

    compare_sampling_strategies(uniform_strategy_cdf = uniform_strategy_cdf,
                                biased_strategy_cdf = biased_strategy_cdf,
                                factors = factors,
                                testing_sample_size = 50 * factors,
                                selection_sample_size = 50 * factors,
                                design_size = factors + 10,
                                regression_formula = model_formulas$linear,
                                filter_thresholds = c(0.1, 0.01, 0.001),
                                gaussian_amplitude = amplitude,
                                gaussian_lengthscale = lengthscale)

    #write.csv(data, "results.csv", row.names = FALSE)
}

test <- run_experiments()
