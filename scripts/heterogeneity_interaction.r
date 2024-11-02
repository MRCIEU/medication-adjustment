library(dplyr)
library(purrr)
library(ggplot2)
library(here)


#' Simulation model
#' @param n sample size
#' @param mean_ldl
#' @param sd_ldl
#' @param b_gy
#' @param b_omed
#' @param med_effect
#' @param abs_adjustment
#' @param rel_adjustment
#' @param sim
#' 
#' @return table of results for estimates on ldl genotype and other genotype under different models
sim <- function(n, mean_ldl=130, sd_ldl=20, b_gy=1, b_omed=1, med_effect=0.8, abs_adjustment = 20, rel_adjustment = 0.8, het_sd = 0, sim=1) {

    args <- c(as.list(environment())) %>% as_tibble()
    g <- rnorm(n)
    g_other <- rnorm(n)

    y <- mean_ldl + g * b_gy + rnorm(n, sd=sqrt(sd_ldl^2-b_gy^2))
    med <- rbinom(n, 1, plogis(scale(y) + g_other * b_omed))

    medeff <- g * med_effect

    y_obs <- y
    y_obs[as.logical(med)] <- y_obs[as.logical(med)] * medeff[as.logical(med)]

    y_adj_abs <- y_obs
    y_adj_abs[as.logical(med)] <- y_obs[as.logical(med)] + abs_adjustment

    y_adj_rel <- y_obs
    y_adj_rel[as.logical(med)] <- y_obs[as.logical(med)] / rel_adjustment


    o <- bind_rows(
        bind_rows(
            summary(lm(y ~ g))$coef[2,],
            summary(lm(y_obs ~ g))$coef[2,],
            summary(lm(y_adj_abs ~ g))$coef[2,],
            summary(lm(y_adj_rel ~ g))$coef[2,],
        ) %>% as_tibble() %>% mutate(measure=c("y", "y_obs", "y_adj_abs", "y_adj_rel"), g="ldl"),

        bind_rows(
            summary(lm(y ~ g_other))$coef[2,],
            summary(lm(y_obs ~ g_other))$coef[2,],
            summary(lm(y_adj_abs ~ g_other))$coef[2,],
            summary(lm(y_adj_rel ~ g_other))$coef[2,],
        ) %>% as_tibble() %>% mutate(measure=c("y", "y_obs", "y_adj_abs", "y_adj_rel"), g="other")
    )
    names(o)[1:4] <- c("beta", "se", "tval", "pval")
    o <- bind_cols(o, args)
    return(o)
}

sim(100000, het_sd=0.2)
