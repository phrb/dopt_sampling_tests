library(Rcpp)
library(rsm)
library(dplyr)
library(AlgDesign)
library(uuid)
library(MASS)
library(logging)
library(stringr)

sourceCpp("rcpp_biased_sample.cpp")
sourceCpp("rcpp_squared_exponential_kernel.cpp")

generate_design_codings <- function(design) {
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

encode_design <- function(design) {
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

biased_sample <- function(sampling_cdf, size, factors, name) {
    rcpp_biased_sample(sampling_cdf,
                       size,
                       factors) %>%
        data.frame() %>%
        mutate(name = name) %>%
        encode_design()
}

filter_formula <- function(complete_formula, p_values, filter_threshold){
    formula_terms <- str_split(paste(deparse(complete_formula), collapse = ""), " ")[[1]]
    filtered_variables <- subset(p_values, values <= filter_threshold)$variables

    filtered_terms <- paste(formula_terms[str_detect(formula_terms,
                                                     paste(filtered_variables,
                                                           collapse = "|"))],
                            collapse = " + ")

    update.formula(complete_formula, paste(". ~ ", filtered_terms, sep = ""))
}

compare_fit <- function(testing_sample, design, formula, filter_thresholds) {
    formula <- formula %>% update.formula(Y ~ .)

    complete_sample <- bind_rows(testing_sample, design)
    anova_summary <- summary(aov(formula, design))

    p_values <- as.data.frame(anova_summary[[1]])["Pr(>F)"]
    names(p_values) <- c("values")

    p_values$variables <- str_replace_all(str_trim(rownames(p_values)),
                                          "I\\(|\\^2\\)",
                                          "")
    rownames(p_values) <- NULL

    results <- NULL

    for(filter_threshold in filter_thresholds) {
        selected <-nrow(subset(p_values, values <= filter_threshold))

        if(is.infinite(filter_threshold) | selected <= 0) {
            loginfo("> Selecting all terms")
            filtered_formula <- formula
        } else {
            loginfo("> Filtering based on p values")
            filtered_formula <- filter_formula(formula, p_values, filter_threshold)
        }

        predictions <- predict(lm(filtered_formula, design), complete_sample)

        fit_distance <- sqrt(sum((complete_sample$Y - predictions) ^ 2))

        min_distance <- abs(unique(min(complete_sample$Y)) -
                            unique(complete_sample$Y[predictions == min(predictions)]))

        if(is.null(results)) {
            results <- data.frame(filter_threshold = filter_threshold,
                                  min_distance = min_distance,
                                  fit_distance = fit_distance)
        } else {
            results <- bind_rows(results, data.frame(filter_threshold = filter_threshold,
                                                     min_distance = min_distance,
                                                     fit_distance = fit_distance))
        }
    }

    return(results)
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
                                         gaussian_lengthscale) {
    loginfo("> Starting to generate samples")

    testing_sample <- biased_sample(uniform_strategy_cdf,
                                    testing_sample_size,
                                    factors,
                                    "testing")

    uniform_selection_sample <- biased_sample(uniform_strategy_cdf,
                                              selection_sample_size,
                                              factors,
                                              "uniform_selection")

    biased_selection_sample <- biased_sample(biased_strategy_cdf,
                                             selection_sample_size,
                                             factors,
                                             "biased_selection")

    loginfo("> Generating design from uniform sample")

    uniform_design <- optFederov(regression_formula,
                                 uniform_selection_sample,
                                 nTrials = design_size)$design %>%
                                                      mutate(name = "uniform")

    loginfo("> Generating design from biased sample")

    biased_design <- optFederov(regression_formula,
                                biased_selection_sample,
                                nTrials = design_size)$design %>%
                                                     mutate(name = "biased")

    loginfo("> Sampling function from gaussian process")

    # Use rbind here to keep encoding information
    sampled_data <- sample_gaussian_process(rbind(testing_sample,
                                                  uniform_design,
                                                  biased_design),
                                            gaussian_amplitude,
                                            gaussian_lengthscale,
                                            -name)

    loginfo("> Comparing uniform fit")

    uniform_fit_comparison <- compare_fit(sampled_data %>%
                                          subset(name == "testing"),
                                          sampled_data %>%
                                          subset(name == "uniform"),
                                          regression_formula,
                                          filter_thresholds)

    uniform_fit_comparison$name <- "uniform"

    loginfo("> Comparing biased fit")

    biased_fit_comparison <- compare_fit(sampled_data %>%
                                         subset(name == "testing"),
                                         sampled_data %>%
                                         subset(name == "biased"),
                                         regression_formula,
                                         filter_thresholds)

    biased_fit_comparison$name <- "biased"

    return(bind_rows(uniform_fit_comparison,
                     biased_fit_comparison))
}

run_experiments <- function(iterations) {
    logging::basicConfig(level = "DEBUG")

    cdf_data <- read.csv("factor_cdfs.csv",
                         header = TRUE)

    factors <- list(linear = 30,
                    quadratic = 30)

    significance_probability <- 0.15
    insignificance_variability <- list(max = 100.0, min = 10.0)
    significance_variability <- list(max = 5.0, min = 0.5)

    amplitude_variability <- list(max = 5.0, min = 1.0)

    strategy_cdfs <- list(uniform = cdf_data %>%
                              subset(name == "uniform") %>%
                              dplyr::select(x, y) %>%
                              as.matrix(),
                          linear = cdf_data %>%
                              subset(name == "linear") %>%
                              dplyr::select(x, y) %>%
                              as.matrix(),
                          quadratic = cdf_data %>%
                              subset(name == "linear") %>%
                              dplyr::select(x, y) %>%
                              as.matrix())

    model_formulas <- list(linear = paste("x", 1:factors$linear, sep = "") %>%
                               paste(sep = "", collapse = " + ") %>%
                               paste("~ ", ., sep = "") %>%
                               formula(),
                           quadratic = paste("x", 1:factors$quadratic, sep = "") %>%
                               paste(sep = "", collapse = " + ") %>%
                               paste("+ ") %>%
                               paste(., paste("I(",
                                              factor_names,
                                              "^2)",
                                              collapse = " + ",
                                              sep = "")) %>%
                               paste("~ ", .) %>%
                               formula())

    design_sizes <- list(linear = factors$linear + 10,
                         quadratic = (2 * factors$quadratic) + 10)

    results <- NULL

    for(model in names(model_formulas)) {
        loginfo(paste("> Generating results for the", model, "model"))

        for(iteration in seq(1, iterations)) {
            loginfo("> Sampling amplitude and lengthscale")

            amplitude <- runif(1,
                               min = amplitude_variability$min,
                               max = amplitude_variability$max)

            lengthscale <- replicate(factors[model][[1]],
                                     ifelse(runif(1) < (1.0 - significance_probability),
                                            runif(1,
                                                  min = insignificance_variability$min,
                                                  max = insignificance_variability$max),
                                            runif(1,
                                                  min = significance_variability$min,
                                                  max = significance_variability$max)))

            loginfo(paste("> Generating results for iteration", iteration))

            result <- compare_sampling_strategies(uniform_strategy_cdf = strategy_cdfs["uniform"][[1]],
                                                  biased_strategy_cdf = strategy_cdfs[model][[1]],
                                                  factors = factors[model][[1]],
                                                  testing_sample_size = 50 * factors[model][[1]],
                                                  selection_sample_size = 50 * factors[model][[1]],
                                                  design_size = design_sizes[model][[1]],
                                                  regression_formula = model_formulas[model][[1]],
                                                  filter_thresholds = c(Inf, 0.1, 0.01, 0.001),
                                                  gaussian_amplitude = amplitude,
                                                  gaussian_lengthscale = lengthscale)

            result$model <- model
            result$iteration <- iteration

            if(is.null(results)) {
                results <- result
            } else {
                results <- bind_rows(results, result)
            }
        }
    }
    #write.csv(data, "results.csv", row.names = FALSE)
    return(results)
}

test <- run_experiments(1)
