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

    medeff <- rnorm(n, med_effect, het_sd)

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



# Setup parameters
params <- expand.grid(
    n = 100000,
    rel_adjustment = seq(0.7, 1, by=0.1),
    het_sd = c(0, 0.1, 0.2),
    sim=1:100
)
dim(params)

# Run analysis on all parameters
res <- purrr::pmap(params, sim, .progress=TRUE) %>% bind_rows()


saveRDS(res, file=here("results/heterogeneity.rds"))

###

res <- readRDS(here("results/heterogeneity.rds"))
res <- res %>% mutate(
    measure = case_when(measure == "y" ~ "True Y", measure == "y_obs" ~ "Observed Y", measure == "y_adj_abs" ~ "Absolute adjustment", measure == "y_adj_rel" ~ "Relative adjustment"),
    g = case_when(g == "ldl" ~ "G->Y", g == "other" ~ "G->Medication")
)


ggplot(res %>% filter(measure != "Absolute adjustment"), aes(x=as.factor(rel_adjustment), y=beta)) +
geom_boxplot(aes(colour=as.factor(het_sd))) +
facet_grid(g ~ measure, scale="free_y") +
labs(y="Genetic effect estimate", x="Relative adjustment to Y", colour="Heterogeneity\nof medication\neffect")
ggsave(here("images/heterogeneity.pdf"))

