library(dplyr)
library(purrr)
library(ggplot2)
library(here)
library(tidyr)

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
sim <- function(n, mean_ldl=130, sd_ldl=20, b_gy=1, b_omed=1, med_effect=0.8, sim=1) {

    args <- c(as.list(environment())) %>% as_tibble()
    g <- rnorm(n)
    g_other <- rnorm(n)

    y <- mean_ldl + g * b_gy + rnorm(n, sd=sqrt(sd_ldl^2-b_gy^2))
    med <- rbinom(n, 1, plogis(scale(y) + g_other * b_omed))
    y_obs <- y
    y_obs[as.logical(med)] <- y_obs[as.logical(med)] * med_effect

    args$b_obs <- lm(y_obs ~ g)$coef[2]
    args$b_cov <- lm(y_obs ~ g + med)$coef[2]
    args$r2_unadj <- cor(y_obs, g)^2

    # range of adjustments
    o1 <- tibble(adjustment=seq(0.6, 1, by=0.01), method="relative", r2=NA, bhat=NA, r_y=NA)
    for(i in 1:nrow(o1)) {
        y_adj_rel <- y_obs
        y_adj_rel[as.logical(med)] <- y_obs[as.logical(med)] / o1$adjustment[i]

        o1$r2[i] <- cor(y_adj_rel, g)^2
        o1$r_y[i] <- cor(y_adj_rel, y)
        o1$bhat[i] <- lm(y_adj_rel ~ g)$coef[2]
    }

    o2 <- tibble(adjustment=seq(10, 40, by=1), method="absolute", r2=NA, bhat=NA, r_y=NA)
    for(i in 1:nrow(o2)) {
        y_adj_abs <- y_obs
        y_adj_abs[as.logical(med)] <- y_obs[as.logical(med)] + o2$adjustment[i]

        o2$r2[i] <- cor(y_adj_abs, g)^2
        o2$r_y[i] <- cor(y_adj_abs, y)
        o2$bhat[i] <- lm(y_adj_abs ~ g)$coef[2]
    }

    o <- bind_cols(bind_rows(o1, o2), args)
    return(o)
}


summarise_sim <- function(o) {
    group_by(o, method) %>%
    arrange(desc(r2)) %>%
    slice_head(n=1)
}


set.seed(12345)
o <- sim(100000, med_effect=0.8)

ggplot(o, aes(x=adjustment, y=r2)) +
geom_line() +
facet_grid(. ~ method, scale="free_x") +
labs(x="Medication adjustment", y="R-square (g, y)")
ggsave(here("images/adjustment_detection_illustration1.pdf"))

ggplot(o, aes(x=adjustment, y=r_y)) +
geom_line() +
facet_grid(. ~ method, scale="free_x")

ggplot(o, aes(y=r2, x=bhat)) +
geom_line() +
facet_grid(. ~ method, scale="free_x") +
geom_vline(xintercept=1, linetype="dotted") +
labs(x="Effect estimate of g on y", y="R-square (g, y)")
ggsave(here("images/adjustment_detection_illustration2.pdf"))


params <- expand.grid(
    n=100000, med_effect=c(0.8),
    sim=1:100
)

fn <- function(...)
{
    sim(...) %>% summarise_sim()
}

res <- pmap(params, fn, .progress=TRUE) %>% bind_rows()
saveRDS(res, file=here("results/adjustment_detection.rds"))

# What fraction of simulations return the correct model?
group_by(res, sim) %>%
arrange(desc(r2)) %>%
slice_head(n=1) %>%
ungroup() %>%
summarise(prel=sum(method=="relative")/n())


res %>% select(-adjustment, -bhat, -r_y, -b_obs) %>% tidyr::pivot_wider(names_from=method, values_from=r2) %>%
ggplot(aes(absolute, relative)) +
geom_point(aes(colour=relative>absolute)) +
geom_abline() +
labs(x="R-square absolute adjustment", y="R-square relative adjustment", colour="Correct model")
ggsave(here("images/adjustment_detection_rsq.pdf"))

group_by(res, sim) %>%
arrange(desc(r2)) %>%
slice_head(n=1) %>%
ungroup() %>% 
mutate(bhat=(1-bhat)^2, b_obs=(1-b_obs)^2) %>%
ggplot(aes(x=bhat, y=b_obs)) +
geom_point(aes(colour=method=="relative")) +
geom_abline() +
labs(x="MSE adjusted", y="MSE unadjusted", colour="Correct model") +
xlim(0,0.25) + ylim(0, 0.25)
ggsave(here("images/adjustment_detection_mse.pdf"))

group_by(res, sim) %>%
arrange(desc(r2)) %>%
slice_head(n=1) %>%
ungroup() %>% 
mutate(bhat=(1-bhat)^2, b_cov=(1-b_cov)^2) %>%
ggplot(aes(x=bhat, y=b_cov)) +
geom_point(aes(colour=method=="relative")) +
geom_abline() +
labs(x="MSE adjusted", y="MSE medication covariate", colour="Correct model") +
xlim(0,0.15) + ylim(0, 0.15)
ggsave(here("images/adjustment_detection_covs_mse.pdf"))

