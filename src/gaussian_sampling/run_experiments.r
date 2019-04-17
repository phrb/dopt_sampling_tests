library(Rcpp)
library(rsm)
library(dplyr)
library(tidyr)
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

filter_formula <- function(complete_formula,
                           p_values,
                           filter_threshold,
                           group_factor_terms = FALSE) {
    formula_terms <- str_split(paste(deparse(complete_formula), collapse = ""), " ")[[1]]

    if(group_factor_terms) {
        p_values$variables <- str_replace_all(p_values$variables,
                                              "I\\(|\\^2\\)|\\^3\\)",
                                              "")
        filtered_variables <- subset(p_values, values <= filter_threshold)$variables

        filtered_terms <- paste(formula_terms[str_detect(formula_terms,
                                                         paste(filtered_variables,
                                                               collapse = "|"))],
                                collapse = " + ")
    } else {
        filtered_terms <- paste(subset(p_values, values <= filter_threshold)$variables,
                                collapse = " + ")
    }

    update.formula(complete_formula, paste(". ~ ", filtered_terms, sep = ""))
}

compute_significance_errors <- function(p_values,
                                        filter_threshold,
                                        gaussian_lengthscale,
                                        gaussian_significance_thresholds) {
    filtered_variables <- unique(subset(p_values, values <= filter_threshold)$variables)

    selected_significant <- data.frame(t(rep(TRUE, length(filtered_variables))))
    names(selected_significant) <- filtered_variables

    real_significant <- data.frame(t(gaussian_lengthscale <=
                                     gaussian_significance_thresholds$min))
    names(real_significant) <- paste("x", 1:length(gaussian_lengthscale), sep = "")

    significances <- bind_rows(real_significant,
                               selected_significant)

    significances[is.na(significances)] <- FALSE

    errors <- data.frame(true_positives = sum(significances[2, ] & significances[1, ]),
                         false_positives = sum(significances[2, ] & !significances[1, ]),
                         false_negatives = sum(!significances[2, ] & significances[1, ]),
                         model_positives = sum(significances[2, ] & TRUE),
                         real_positives = sum(significances[1, ] & TRUE))
}

compare_fit <- function(testing_sample,
                        design,
                        formula,
                        filter_thresholds,
                        noise_sd,
                        gaussian_lengthscale,
                        gaussian_significance_thresholds) {
    formula <- formula %>% update.formula(Y ~ .)

    design_with_noise <- design
    design_with_noise$Y <- design_with_noise$Y + rnorm(1, mean = 0.0, sd = noise_sd)

    anova_summary <- summary(aov(formula, design_with_noise))

    complete_sample <- bind_rows(testing_sample, design)

    p_values <- as.data.frame(anova_summary[[1]])["Pr(>F)"]
    names(p_values) <- c("values")

    p_values$variables <- str_trim(rownames(p_values))
    rownames(p_values) <- NULL

    lengthscale_data <- data.frame(t(gaussian_lengthscale))
    names(lengthscale_data) <- paste("lengthscale_x", 1:length(gaussian_lengthscale), sep = "")

    results <- NULL


    for(filter_threshold in filter_thresholds) {
        selected <-nrow(subset(p_values, values <= filter_threshold))

        if(is.infinite(filter_threshold) | selected <= 0) {
            filtered_formula <- formula
        } else {
            filtered_formula <- filter_formula(formula, p_values, filter_threshold)
        }

        predictions <- predict(lm(filtered_formula, design_with_noise), complete_sample)

        fit_distance <- sqrt(sum((complete_sample$Y - predictions) ^ 2))

        min_distance <- abs(min(complete_sample$Y) -
                            complete_sample$Y[predictions == min(predictions)][1])


        errors <- compute_significance_errors(p_values,
                                              filter_threshold,
                                              gaussian_lengthscale,
                                              gaussian_significance_thresholds)

        result <- data.frame(filter_threshold = filter_threshold,
                             min_distance = min_distance,
                             fit_distance = fit_distance,
                             true_positives = errors$true_positives,
                             false_positives = errors$false_positives,
                             false_negatives = errors$false_negatives,
                             model_positives = errors$model_positives,
                             real_positives = errors$real_positives)

        result <- bind_cols(result, lengthscale_data)

        if(is.null(results)) {
            results <- result
        } else {
            results <- bind_rows(results, result)
        }
    }

    return(results)
}

compare_sampling_strategies <- function(uniform_strategy_cdf,
                                        biased_strategy_cdf,
                                        factors,
                                        testing_sample_size,
                                        selection_sample_size,
                                        design_size,
                                        regression_formula,
                                        filter_thresholds,
                                        noise_sd,
                                        gaussian_amplitude,
                                        gaussian_lengthscale,
                                        gaussian_significance_thresholds) {
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
                                          filter_thresholds,
                                          noise_sd,
                                          gaussian_lengthscale,
                                          gaussian_significance_thresholds)

    uniform_fit_comparison$name <- "uniform"

    loginfo("> Comparing biased fit")

    biased_fit_comparison <- compare_fit(sampled_data %>%
                                         subset(name == "testing"),
                                         sampled_data %>%
                                         subset(name == "biased"),
                                         regression_formula,
                                         filter_thresholds,
                                         noise_sd,
                                         gaussian_lengthscale,
                                         gaussian_significance_thresholds)

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

    model_size <- list(linear = 30,
                       quadratic = 60)

    filter_thresholds <- c(Inf, 0.1, 0.01, 0.001)
    noise_sd <- 4.0

    significance_probability <- 0.1
    insignificance_variability <- list(max = 10.0, min = 5.0)
    significance_variability <- list(max = 1.1, min = 0.9)

    significance_thresholds <- list(min = significance_variability$max,
                                    max = insignificance_variability$min)

    amplitude_variability <- list(max = 1.5, min = 1.0)

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
                               paste(., paste("I(x",
                                              1:factors$quadratic,
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
                                                  testing_sample_size = 100 * factors[model][[1]],
                                                  selection_sample_size = 100 * factors[model][[1]],
                                                  design_size = design_sizes[model][[1]],
                                                  regression_formula = model_formulas[model][[1]],
                                                  filter_thresholds = filter_thresholds,
                                                  noise_sd = noise_sd,
                                                  gaussian_amplitude = amplitude,
                                                  gaussian_lengthscale = lengthscale,
                                                  gaussian_significance_thresholds = significance_thresholds)

            result$model <- model
            result$id <- UUIDgenerate()
            result$factors <- factors[model][[1]]
            result$model_size <- model_size[model][[1]]
            result$amplitude <- amplitude

            if(is.null(results)) {
                results <- result
            } else {
                results <- bind_rows(results, result)
            }
        }
    }

    write.csv(results, "results.csv", row.names = FALSE)
    return(results)
}

# test <- run_experiments(1)
run_experiments(10)
