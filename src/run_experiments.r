library(rsm)
library(dplyr)
library(AlgDesign)
library(uuid)

factor_cdf <- function(f, name = "Generic Factors", interval = 0.0001) {
    x <- seq(0.0, 1.0, interval)
    y <- mapply(f, x)

    data <- data.frame(x = x, y = y, name = name)
    factor_pdf <- data.frame(x = x, y = data$y, name = paste(name, "PDF"))

    factor_pdf$y <- factor_pdf$y / sum(factor_pdf$y)
    data <- bind_rows(data, factor_pdf)

    y_cdf <- c()

    target_name <- paste(name, "PDF")

    pdf_data <- subset(data, name == target_name)$y

    for (i in 1:length(pdf_data)) {
        y_cdf <- c(y_cdf, sum(pdf_data[1:i]))
    }

    target_name <- paste(name, "CDF")
    data <- bind_rows(data,
                      data.frame(x = x, y = y_cdf, name = target_name))
    return(data)
}

sample_factor <- function(cdf_data, samples) {
    x_sampled <- c()

    for (i in 1:samples) {
        new_sample <- cdf_data$x[sum(cdf_data$y <= runif(1)) + 1]
        x_sampled <- c(x_sampled, new_sample)
    }

    return(x_sampled)
}

get_runif_design_codings <- function(design) {
    design_codings <- list()

    for (i in 1:ncol(design)) {
        design_codings[[paste("x", i, sep = "")]] <- formula(paste(paste("x", i, "~", sep = ""),
                                                                  paste("(X", i, "-0.5)/0.5", sep = ""),
                                                                  sep = ""))
    }

    return(design_codings)
}

get_runif_coded_design <- function(design) {
    design_codings <- get_runif_design_codings(design)
    return(coded.data(design, formulas = design_codings))
}

cubic_function <- function(x) {
    return((64 * (0.5 - x) ** 6) + (dnorm(x, 0.35, 0.04) / 10) + (dnorm(x, 0.65, 0.04) / 10))
}

quadratic_function <- function(x) {
    return((64 * (0.5 - x) ** 6) + (dnorm(x, 0.5, 0.04) / 10))
}

linear_function <- function(x) {
    return((64 * (0.5 - x) ** 6))
}

compute_response <- function(coefficients, formula, experiment) {
    return(as.vector((model.matrix(formula, experiment) %*% coefficients)))
}

add_noise <- function(design, sd = 1) {
    return(as.vector(design$Y + rnorm(nrow(design), mean = 0, sd = sd)))
}

compare_fit <- function(coefficients, sample, formula, response_formula, design_a, design_b, sample_size = 100, sd = 1) {
    return(c(sqrt(sum(((compute_response(coefficients, formula, sample)) -
                        predict(lm(response_formula, design_a),
                               sample)) ^ 2)),
             sqrt(sum(((compute_response(coefficients, formula, sample)) -
                        predict(lm(response_formula, design_b),
                               sample)) ^ 2))))
}

compare_post_anova_fit <- function(coefficients,
                                   sample,
                                   formula,
                                   response_formula,
                                   design_a,
                                   design_b,
                                   sample_size = 100,
                                   sd = 1,
                                   prf_threshold = 0.0001) {
    regression_summary_a <- summary(aov(response_formula, design_a))

    prf_values_a <- as.data.frame(regression_summary_a[[1]])["Pr(>F)"]
    names(prf_values_a) <- c("values")

    if(nrow(subset(prf_values_a, values <= prf_threshold)) > 0) {
        refitted_formula_a <- formula(paste("Y ~", paste(rownames(subset(prf_values_a, values <= prf_threshold)), sep = "", collapse = " + ")))
    } else {
        refitted_formula_a <- response_formula
    }

    regression_summary_b <- summary(aov(response_formula, design_b))

    prf_values_b <- as.data.frame(regression_summary_b[[1]])["Pr(>F)"]
    names(prf_values_b) <- c("values")

    if(nrow(subset(prf_values_b, values <= prf_threshold)) > 0) {
      refitted_formula_b <- formula(paste("Y ~", paste(rownames(subset(prf_values_b, values <= prf_threshold)), sep = "", collapse = " + ")))
    } else {
      refitted_formula_b <- response_formula
    }

    return(c(sqrt(sum(((compute_response(coefficients, formula, sample)) -
                       predict(lm(refitted_formula_a, design_a),
                               sample)) ^ 2)),
            sqrt(sum(((compute_response(coefficients, formula, sample)) -
                      predict(lm(refitted_formula_b, design_b),
                              sample)) ^ 2))))
}

# https://stats.stackexchange.com/questions/93540/testing-equality-of-coefficients-from-two-different-regressions
#return(sum((coefficients - unname(coef(lm(Y ~., design)))) / sqrt((data.frame(summary(lm(Y ~., design))["coefficients"])[, 2] ^ 2))))
compare_coefficients <- function(coefficients, formula, design) {
    return(sqrt(sum((coefficients - unname(coef(lm(formula, design)))) ^ 2)))
}

compare_post_anova_coefficients <- function(coefficients,
                                            formula,
                                            design_a,
                                            design_b,
                                            prf_threshold = 0.0001) {
    regression_summary_a <- summary(aov(formula, design_a))

    prf_values_a <- as.data.frame(regression_summary_a[[1]])["Pr(>F)"]
    names(prf_values_a) <- c("values")

    if(nrow(subset(prf_values_a, values <= prf_threshold)) > 0) {
      refitted_formula_a <- formula(paste("Y ~", paste(rownames(subset(prf_values_a, values <= prf_threshold)), sep = "", collapse = " + ")))
    } else {
      refitted_formula_a <- formula
    }

    regression_summary_b <- summary(aov(formula, design_b))

    prf_values_b <- as.data.frame(regression_summary_b[[1]])["Pr(>F)"]
    names(prf_values_b) <- c("values")

    if(nrow(subset(prf_values_b, values <= prf_threshold)) > 0) {
      refitted_formula_b <- formula(paste("Y ~", paste(rownames(subset(prf_values_b, values <= prf_threshold)), sep = "", collapse = " + ")))
    } else {
      refitted_formula_b <- formula
    }

    named_coefficients <- data.frame(t(coefficients))
    names(named_coefficients) <- names(coef(lm(formula, design_a)))

    named_refitted_coefficients_a <- data.frame(t(coef(lm(refitted_formula_a, design_a))))
    names(named_refitted_coefficients_a) <- names(coef(lm(refitted_formula_a, design_a)))

    named_refitted_coefficients_b <- data.frame(t(coef(lm(refitted_formula_b, design_b))))
    names(named_refitted_coefficients_b) <- names(coef(lm(refitted_formula_b, design_b)))

    coefficient_data <- bind_rows(named_coefficients,
                                  named_refitted_coefficients_a,
                                  named_refitted_coefficients_b)

    coefficient_data[is.na(coefficient_data)] <- 0.0

    return(c(sqrt(sum((as.numeric(coefficient_data[1, ]) - as.numeric(coefficient_data[2, ])) ^ 2)),
             sqrt(sum((as.numeric(coefficient_data[1, ]) - as.numeric(coefficient_data[3, ])) ^ 2))))
}

compare_real_coefficients <- function(coefficients, coefficient_limits, formula, design) {
    prf_values <- unname(summary(lm(formula, design))[4][["coefficients"]][,4])
    return(c(sum(prf_values[abs(coefficients) > coefficient_limits$min])))
}

compare_noise_coefficients <- function(coefficients, coefficient_limits, formula, design) {
    prf_values <- unname(summary(lm(formula, design))[4][["coefficients"]][,4])
    return(sum(prf_values[abs(coefficients) <= coefficient_limits$min]))
}

linear_experiment <- function(coefficient_number,
                              coefficient_probability,
                              coefficient_variability,
                              coefficient_noise,
                              factors,
                              fit_comparison_sample_size,
                              federov_samples,
                              formula,
                              response_formula,
                              samples,
                              lin_cdf,
                              noise_sd_vector,
                              prf_threshold_vector) {

    experiments <- NULL
    coefficients <- replicate(coefficient_number,
                              ifelse(runif(1) > (1.0 - coefficient_probability),
                              ifelse(runif(1) >= 0.5,
                                     runif(1,
                                           min = -coefficient_variability$max,
                                           max = -coefficient_variability$min),
                                     runif(1,
                                           min = coefficient_variability$min,
                                           max = coefficient_variability$max)),
                              rnorm(1, mean = 0, sd = coefficient_noise)))

    current_uuid = UUIDgenerate()

    for (prf_threshold in prf_threshold_vector) {
        for (noise_sd in noise_sd_vector) {
            fit_comparison_sample <- get_runif_coded_design(data.frame(replicate(factors,
                                                                                 runif(fit_comparison_sample_size))))


            federov_sample <- get_runif_coded_design(data.frame(replicate(factors,
                                                                          runif(federov_samples))))

            federov_design <- optFederov(formula,
                                         federov_sample,
                                         nTrials = samples)$design

            federov_design$Y <- compute_response(coefficients, formula, federov_design)
            federov_design$Y <- add_noise(federov_design, sd = noise_sd)

            federov_sample <- get_runif_coded_design(data.frame(replicate(factors,
                                                                          sample_factor(lin_cdf,
                                                                                        federov_samples))))

            biased_federov_design <- optFederov(formula,
                                                federov_sample,
                                                nTrials = samples)$design

            biased_federov_design$Y <- compute_response(coefficients, formula, biased_federov_design)
            biased_federov_design$Y <- add_noise(biased_federov_design, sd = noise_sd)

            coefficient_distance = c(compare_coefficients(coefficients,
                                                          response_formula,
                                                          federov_design),
                                     compare_coefficients(coefficients,
                                                          response_formula,
                                                          biased_federov_design))

            coefficient_refit_distance = compare_post_anova_coefficients(coefficients,
                                                                         response_formula,
                                                                         federov_design,
                                                                         biased_federov_design,
                                                                         prf_threshold = prf_threshold)

            fit_distance = compare_fit(coefficients,
                                       fit_comparison_sample,
                                       formula,
                                       response_formula,
                                       federov_design,
                                       biased_federov_design,
                                       sample_size = fit_comparison_sample_size,
                                       sd = noise_sd)

            refit_distance = compare_post_anova_fit(coefficients,
                                                    fit_comparison_sample,
                                                    formula,
                                                    response_formula,
                                                    federov_design,
                                                    biased_federov_design,
                                                    sample_size = fit_comparison_sample_size,
                                                    sd = noise_sd,
                                                    prf_threshold = prf_threshold)

            new_experiment <- data.frame(D                          = c(eval.design(formula, federov_design)$determinant,
                                                                        eval.design(formula, biased_federov_design)$determinant),
                                         coefficient_distance       = coefficient_distance,
                                         coefficient_refit_distance = coefficient_refit_distance,
                                         fit_distance               = fit_distance,
                                         refit_distance             = refit_distance,
                                         names                      = c("Federov with Uniform Sample", "Federov with Biased Sample"),
                                         id                         = current_uuid,
                                         prft                       = prf_threshold,
                                         noise                      = noise_sd,
                                         coefficient_probability    = coefficient_probability,
                                         regression_model           = "linear",
                                         normal_samples             = 0.0,
                                         coefficients               = coefficient_number)

            if (is.null(experiments)) {
                experiments <- new_experiment
            } else {
                experiments <- bind_rows(experiments, new_experiment)
            }
        }
    }

    return(experiments)
}

quadratic_experiment <- function(coefficient_number,
                                 coefficient_probability,
                                 coefficient_variability,
                                 coefficient_noise,
                                 factors,
                                 fit_comparison_sample_size,
                                 federov_samples,
                                 formula,
                                 response_formula,
                                 samples,
                                 quad_cdf,
                                 normal_samples_ratio,
                                 noise_sd_vector,
                                 prf_threshold_vector) {

    experiments <- NULL
    coefficients <- replicate(coefficient_number,
                              ifelse(runif(1) > (1.0 - coefficient_probability),
                              ifelse(runif(1) >= 0.5,
                                     runif(1,
                                           min = -coefficient_variability$max,
                                           max = -coefficient_variability$min),
                                     runif(1,
                                           min = coefficient_variability$min,
                                           max = coefficient_variability$max)),
                              rnorm(1, mean = 0, sd = coefficient_noise)))

    current_uuid = UUIDgenerate()

    for (prf_threshold in prf_threshold_vector) {
        for (noise_sd in noise_sd_vector) {
            fit_comparison_sample <- get_runif_coded_design(data.frame(replicate(factors,
                                                                                 runif(fit_comparison_sample_size))))


            federov_sample <- get_runif_coded_design(data.frame(replicate(factors,
                                                                          runif(federov_samples))))

            federov_design <- optFederov(formula,
                                         federov_sample,
                                         nTrials = samples)$design

            federov_design$Y <- compute_response(coefficients, formula, federov_design)
            federov_design$Y <- add_noise(federov_design, sd = noise_sd)

            federov_sample <- get_runif_coded_design(data.frame(replicate(factors,
                                                                          sample_factor(quad_cdf,
                                                                                        as.integer((1 - normal_samples_ratio) * federov_samples)))))

            if(normal_samples_ratio > 0) {
                federov_sample <- rbind(federov_sample, get_runif_coded_design(data.frame(replicate(factors,
                                                                                                    rnorm(as.integer(normal_samples_ratio * federov_samples),
                                                                                                          mean = 0.5, sd = 0.04)))))
            }

            biased_federov_design <- optFederov(formula,
                                                federov_sample,
                                                nTrials = samples)$design

            biased_federov_design$Y <- compute_response(coefficients, formula, biased_federov_design)
            biased_federov_design$Y <- add_noise(biased_federov_design, sd = noise_sd)

            coefficient_distance = c(compare_coefficients(coefficients,
                                                          response_formula,
                                                          federov_design),
                                     compare_coefficients(coefficients,
                                                          response_formula,
                                                          biased_federov_design))

            coefficient_refit_distance = compare_post_anova_coefficients(coefficients,
                                                                         response_formula,
                                                                         federov_design,
                                                                         biased_federov_design,
                                                                         prf_threshold = prf_threshold)

            fit_distance = compare_fit(coefficients,
                                       fit_comparison_sample,
                                       formula,
                                       response_formula,
                                       federov_design,
                                       biased_federov_design,
                                       sample_size = fit_comparison_sample_size,
                                       sd = noise_sd)

            refit_distance = compare_post_anova_fit(coefficients,
                                                    fit_comparison_sample,
                                                    formula,
                                                    response_formula,
                                                    federov_design,
                                                    biased_federov_design,
                                                    sample_size = fit_comparison_sample_size,
                                                    sd = noise_sd,
                                                    prf_threshold = prf_threshold)

            new_experiment <- data.frame(D                          = c(eval.design(formula, federov_design)$determinant,
                                                                        eval.design(formula, biased_federov_design)$determinant),
                                         coefficient_distance       = coefficient_distance,
                                         coefficient_refit_distance = coefficient_refit_distance,
                                         fit_distance               = fit_distance,
                                         refit_distance             = refit_distance,
                                         names                      = c("Federov with Uniform Sample", "Federov with Biased Sample"),
                                         id                         = current_uuid,
                                         prft                       = prf_threshold,
                                         noise                      = noise_sd,
                                         coefficient_probability    = coefficient_probability,
                                         regression_model           = "quadratic",
                                         normal_samples             = normal_samples_ratio,
                                         coefficients               = coefficient_number)

            if (is.null(experiments)) {
                experiments <- new_experiment
            } else {
                experiments <- bind_rows(experiments, new_experiment)
            }
        }
    }

    return(experiments)
}

run_experiments <- function() {
    lin_data <- factor_cdf(linear_function, name = "Linear Factors")
    lin_cdf <- subset(lin_data, name == "Linear Factors CDF")

    quad_data <- factor_cdf(quadratic_function, name = "Quadratic Factors")
    quad_cdf <- subset(quad_data, name == "Quadratic Factors CDF")

    data <- NULL

    #runs <- 100
    runs <- 1
    factors <- 60

    federov_samples <- 1000
    fit_comparison_sample_size <- 1000

    normal_samples_sd <- 0.04
    #noise_sd <- c(1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
    #prf_threshold <- c(0.1, 0.01, 0.001)

    noise_sd <- c(4.0)
    prf_threshold <- c(0.01)

    coefficient_probability <- 0.15
    coefficient_variability <- list(max = 7, min = 1)
    coefficient_noise = 0.01

    factor_names <- paste("x", 1:factors, sep = "")

    samples <- list(linear = (1 * factors) + 10,
                    quadratic = (2 * factors) + 10)

    coefficient_number <- list(linear = factors + 1,
                               quadratic = (2 * factors) + 1)

    normal_samples_ratio <- list(linear = 0.08,
                                 quadratic = 0.0)

    model_formula <- list(linear = formula(paste("~ ",
                                                 paste(factor_names,
                                                       sep = "",
                                                       collapse = " + "),
                                                 sep = "")),
                          quadratic = formula(paste("~ ",
                                                    paste(factor_names,
                                                          sep = "",
                                                          collapse = " + "),
                                                    " + ",
                                                    paste("I(",
                                                          factor_names,
                                                          "^2)",
                                                          sep = "",
                                                          collapse = " + "),
                                                    sep = "")))

    response_formula <- list(linear = formula(paste("Y ~ ",
                                                    paste(factor_names,
                                                          sep = "",
                                                          collapse = " + "),
                                                    sep = "")),
                             quadratic = formula(paste("Y ~ ",
                                                       paste(factor_names,
                                                             sep = "",
                                                             collapse = " + "),
                                                       " + ",
                                                       paste("I(",
                                                             factor_names,
                                                             "^2)",
                                                             sep = "",
                                                             collapse = " + "),
                                                       sep = "")))

    for (i in 1:runs) {
        if (i %% 10 == 0) {
            print(i)
        }

        new_data <- linear_experiment(coefficient_number$linear,
                                      coefficient_probability,
                                      coefficient_variability,
                                      coefficient_noise,
                                      factors,
                                      fit_comparison_sample_size,
                                      federov_samples,
                                      model_formula$linear,
                                      response_formula$linear,
                                      samples$linear,
                                      lin_cdf,
                                      noise_sd,
                                      prf_threshold)

        new_data <- bind_rows(new_data, quadratic_experiment(coefficient_number$quadratic,
                                                             coefficient_probability,
                                                             coefficient_variability,
                                                             coefficient_noise,
                                                             factors,
                                                             fit_comparison_sample_size,
                                                             federov_samples,
                                                             model_formula$quadratic,
                                                             response_formula$quadratic,
                                                             samples$quadratic,
                                                             quad_cdf,
                                                             normal_samples_ratio$quadratic,
                                                             noise_sd,
                                                             prf_threshold))

        if (is.null(data)) {
            data <- new_data
        } else {
            data <- bind_rows(data, new_data)
        }
    }

    return(data)
}

data <- run_experiments()
write.csv(data, "results.csv", row.names = FALSE)
