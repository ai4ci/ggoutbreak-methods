# Supplementary material.
# install.packages("ggtext")

library(tidyverse)
library(ggoutbreak)

.gg_pedantic()

# We are interested in the probability integral transform histogram.
# suppose we have a set of true values.

rlnorm4 = function(n, median, coef_variation) {
  mu = log(median)
  sigma = sqrt(log(coef_variation^2 + 1))
  return(rlnorm(n, mu, sigma))
}

true_value = rnorm(1000, 0, 1)
n = length(true_value)

# suppose we have an set of probabilistic estimates for the true values,
# based on some algorithm which has an error associated with it:
# Because the error function has a non zero mean the overall estimate is biased

est_mean = true_value + rnorm(1000, mean = 0.2, sd = 1)

# The SD of our estimate (related to the confidence intervals) is compensating
# for the sd of our error function. In this case it is smaller than the sd of
# our error function.
est_sd = rlnorm4(1000, median = 0.5, coef_variation = 0.1)
# and lets assume our probabilistic estimate is normally distributed

# the quantiles of the true value as judged by our estimates:
pit_value = pnorm(true_value, est_mean, est_sd)

# This is the PIT histogram which is U shaped:
# This is because our estimates noise SD (1) is larger than our estimates average SD (0.25)
# => our predictions are over confident and the true value is outside the ends
# of our predictions most of the time.
ggplot(tibble(pit = pit_value), aes(x = pit)) + geom_histogram(binwidth = 0.01)

# We are looking for a metric to quantify the level of inaccuracy numerically.
# The ideal for a PIT histogram is the uniform. So we think of a wasserstein
# distance.

# if each pit value (sorted ascending) is moved to its correct place in a
# uniform line and we take the sum of the absolute values of this movement gives
# the wasserstein distance.

# if we had 10 items the correct place though would be seq(0.05,0.95,length.out=10)

unif = seq(1 / (2 * n), 1 - 1 / (2 * n), length.out = n)
was = sum(abs(unif - sort(pit_value)))
# we can normalise this by the number of estimates:
was = was / n

# if there is more mass being moved to the left, compared to moved to the right
# we have an indicator of degree of bias.
was_bias = sum(unif - sort(pit_value))
was_bias = was_bias / n

# However given we have estimate means and true values we have an direct estimate of
# the bias
bias = mean(est_mean - true_value)

# If we unbias out estimates we can recalculate pit
unbiased_pit_value = pnorm(true_value, est_mean - bias, est_sd)
unbiased_was = mean(abs(unif - sort(unbiased_pit_value)))

ggplot(tibble(pit = unbiased_pit_value), aes(x = pit)) +
  geom_histogram(binwidth = 0.01)


# We really want to also include in our metric distance from the towards or
# away from the centre. A positive value here is towards the centre and
# represents overconfidence.
directed_was = mean(
  (unif - sort(unbiased_pit_value)) * ifelse(unif < 0.5, 1, -1)
)

# calculate a CRPS for one probabilistic estimate
crps = function(true_value, cdf, min, max) {
  x = true_value
  fn = rlang::as_function(cdf)
  tryCatch(
    stats::integrate(
      f = \(y) fn(y)^2,
      lower = min,
      upper = x,
      stop.on.error = FALSE
    )$value +
      +stats::integrate(
        f = \(y) (fn(y) - 1)^2,
        lower = x,
        upper = max,
        stop.on.error = FALSE
      )$value,
    error = function(e) {
      NA
    }
  )
}

crps_norm = function(true_value, mean, sd) {
  z = (true_value - mean) / sd
  sd * (2 * dnorm(z) + z * (2 * pnorm(z) - 1) - 2 / sqrt(pi))
}

crps_norm(0, 0, 1)
scoringRules::crps_cnorm(0, 0, 1)

# build this into a purrr pipeline
# each row here represents one estimator with error function and confidence
# parameters
# each estimator is assessed against one list of (at the moment 1000) true values
# The ground truth distribution in N(0,1) but this is arbitrary as only defines
# the offset of the estimation problem, all the analysis here is based on error
# and bias of the estimator. We could have used 0 for all the true values.
#
#
run_scenario = function(true_value, param_df, estfn = pnorm) {
  # param_df has columns est_error_mean, est_error_sd, conf_sd_median, conf_sd_cv = 0,
  n = length(true_value)
  unif = seq(1 / (2 * n), 1 - 1 / (2 * n), length.out = n)

  if (!"conf_sd_cv" %in% colnames(param_df)) {
    param_df = param_df %>% mutate(conf_sd_cv = 0)
  }

  # error function is assumed to be normally distributed, with mean = bias, sd = sd / noise.
  # confidence distribution will have a
  param_df = param_df %>%
    mutate(
      est_mean = purrr::map2(
        est_error_mean,
        est_error_sd,
        ~ true_value + rnorm(n, mean = .x, sd = .y)
      ),
      est_sd = purrr::map2(
        conf_sd_median,
        conf_sd_cv,
        ~ rlnorm4(n, median = .x, coef_variation = .y)
      ),
    )

  param_df = param_df %>%
    mutate(
      # the quantiles of the true value as judged by our estimates:
      pit_value = purrr::map2(est_mean, est_sd, ~ estfn(true_value, .x, .y)),
      percent_95_ci = purrr::map_dbl(
        pit_value,
        ~ mean(.x > 0.025 & .x < 0.975)
      ),
      percent_50_ci = purrr::map_dbl(pit_value, ~ mean(.x > 0.25 & .x < 0.75)),
      crps = purrr::map2(
        est_mean,
        est_sd,
        ~ scoringRules::crps_norm(true_values, .x, .y)
      ),
      mean_crps = purrr::map_dbl(crps, mean),
      # The absolute bias in the estimates
      # Could this have been done on median rather than mean. Median may be more often
      # available (could solve CDF for 0.5) rather than integrate PDF (which we don;t have) for mean.
      # There is rationale to use median. and it happens to be the same for the
      # normal distribution.
      bias = purrr::map_dbl(est_mean, ~ mean(.x - true_value)),
      unbiased = purrr::map2(est_mean, bias, ~ .x - .y), #est mean is a column of vectors
      unbiased_pit_value = purrr::map2(
        unbiased,
        est_sd,
        ~ estfn(true_value, .x, .y)
      ),
      unbiased_percent_95_ci = purrr::map_dbl(
        unbiased_pit_value,
        ~ mean(.x > 0.025 & .x < 0.975)
      ),
      unbiased_percent_50_ci = purrr::map_dbl(
        unbiased_pit_value,
        ~ mean(.x > 0.25 & .x < 0.75)
      ),
      unbiased_crps = purrr::map2(
        unbiased,
        est_sd,
        ~ scoringRules::crps_norm(true_values, .x, .y)
      ),
      unbiased_mean_crps = purrr::map_dbl(unbiased_crps, mean),
      # the wasserstein distance.
      was = purrr::map_dbl(pit_value, ~ mean(abs(unif - sort(.x)))),
      # the bias.
      was_bias = purrr::map_dbl(pit_value, ~ mean(unif - sort(.x))),
      # the bias corrected wasserstein distance.
      unbiased_was = purrr::map_dbl(
        unbiased_pit_value,
        ~ mean(abs(unif - sort(.x)))
      ),
      # the bias corrected wasserstein distance with direction from centre.
      directed_was = purrr::map_dbl(
        unbiased_pit_value,
        ~ mean((unif - sort(.x)) * ifelse(unif < 0.5, 1, -1))
      ),
      label = sprintf(
        "PIT-W: %1.3f\nadj PIT-W: %1.3f\ndir PIT-W: %1.3f\nmean CRPS: %1.3f",
        was,
        unbiased_was,
        directed_was,
        mean_crps
      ),
      settings = sprintf(
        "Err=N(%1.2g, %1.2g); Conf SD=%1.2g",
        est_error_mean,
        est_error_sd,
        conf_sd_median
      ),
      confidence = case_when(
        est_error_sd > conf_sd_median ~ "overconfident",
        est_error_sd < conf_sd_median ~ "conservative",
        TRUE ~ "appropriate"
      ) %>%
        factor(levels = c("conservative", "appropriate", "overconfident")),

      biased = case_when(
        est_error_mean > 0 ~ "positively",
        est_error_mean < 0 ~ "negatively",
        TRUE ~ "unbiased"
      ) %>%
        factor(levels = c("negatively", "unbiased", "positively"))
    )

  return(param_df)
}


# Figure 2 ----

true_values = rnorm(10000, 0, 1)
# tmp = run_scenario( true_values, tibble(est_error_mean = 0.5, est_error_sd = 1.5, conf_sd_median = 1) )

plot_scenarios = tidyr::crossing(
  est_error_mean = c(-0.5, 0, 0.5),
  est_error_sd = c(1 / 1.5, 1, 1.5),
  conf_sd_median = 1
)

plot_data = run_scenario(true_values, plot_scenarios)
hist_data = plot_data %>%
  select(confidence, biased, pit_value) %>%
  unnest(c(pit_value))
fig1 = ggplot() +
  geom_histogram(
    data = hist_data,
    mapping = aes(x = pit_value),
    binwidth = 0.02
  ) +
  ggpp::geom_label_npc(
    data = plot_data,
    aes(label = label),
    npcy = 0.95,
    npcx = 0.05,
    vjust = "top",
    hjust = "left",
    size = 4,
    size.unit = "pt",
    alpha = 0.75,
    label.r = unit(0, "pt")
  ) +
  facet_grid(confidence ~ biased, scales = "free_y") +
  xlab("PIT value")

.gg_save_as(
  fig1,
  here::here("s2/fig/fig2-pit-histograms"),
  size = std_size$half,
  formats = c("pdf", "eps")
)

plot_data %>%
  transmute(
    confidence,
    biased,
    `bias (µₑ)` = sprintf("%1.2f", est_error_mean),
    `sharpness (σₑ)` = sprintf("%1.2f", conf_sd_median),
    `calibration (σₑ/σₛ)` = sprintf("%1.2f", conf_sd_median / est_error_sd),
    `kappa (κ)` = 0
  ) %>%
  .hux_tidy(rowGroupVars = c("biased", "confidence"), colGroupVars = c()) %>%
  .hux_save_as("s2/fig/tab1-params.pdf", size = std_size$third)

# Doesn't much matter if you fix the estimate error and vary the confidence
# plot_scenarios = tidyr::crossing(
#   est_error_mean = c(-0.5,0,0.5),
#   conf_sd_mean = c(1/1.5,1,1.5),
#   est_error_sd = 1
# )

# Figure 3 ----

bias_limit = 1
conf_range = 2

analysis_scenarios = tibble(
  sd_sd_mean = c(conf_range),
  sd_label = forcats::as_factor(c("heterogenous CIs"))
) %>%
  group_by(sd_label, sd_sd_mean) %>%
  reframe(
    est_error_mean = runif(1000, -bias_limit, bias_limit),
    est_error_sd = exp(runif(1000, -log(conf_range), log(conf_range))),
    conf_sd_median = exp(runif(1000, -log(conf_range), log(conf_range))),
    conf_sd_cv = runif(1000, 0, sd_sd_mean)
  ) %>%
  ungroup()
analysis_data = run_scenario(true_value, analysis_scenarios) %>%
  mutate(
    # convert mean to median assuming uninformed beta distribution prior to
    # allow values with zero counts to be represented
    unbiased_percent_50_ci = (unbiased_percent_50_ci * n + 1 / 3) / (n + 2 / 3),
    conf_vs_error = conf_sd_median / est_error_sd
  )

# estimate SD/error SD:
# ratio between the SD of an estimate and the SD of the error function
# a low value means that the estimate is more confident than it should be

# directed PIT wasserstein
# a bias corrected measure.
# designed to show PIT EMD distance from a uniform distribution with direction
# towards or away from the centre

# bias adjusted coverage probability
# if the probabilistic estimate is adjusted to correct for average bias for that
# estimator, would the true value have fallen in the 95% CI for the estimator?
# replicated over multiple estimations what is the probability.
# If the estimate is narrow this will be low. If the estimate is conservative this
# will be high (>0.95)
# If the variability of the width of the estimator is high for some reason the
# coverage probability only goes down. This is because if there is a spread of
# CIs, the more uncertain CIs will catch everything fixed ones can, but more
# certain CIs will more likely result in no coverage. Because it is binary the
# distribution of mass behaviour is complex

p2a = ggplot(
  analysis_data,
  aes(x = conf_vs_error, y = mean_crps, colour = conf_sd_cv)
) +
  geom_point() +
  scale_x_log10() +
  # facet_wrap(~sd_label)+
  scale_color_viridis_c(name = "CI variation") +
  geom_vline(xintercept = 1) +
  # ggpp::annotate(geom=ggpp::GeomLabelNpc, label="conservative \u2192", npcx=0.6,npcy=0.95,hjust="left",vjust="center",alpha=0.75, size=6, size.unit="pt", label.r=unit(0,"pt"))+
  # ggpp::annotate(geom=ggpp::GeomLabelNpc, label="\u2190 overconfident", npcx=0.4,npcy=0.95,hjust="right",vjust="center",alpha=0.75, size=6, size.unit="pt", label.r=unit(0,"pt"))+
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = Inf,
    hjust = -0.5,
    vjust = 1.5,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = Inf,
    hjust = 1.5,
    vjust = 1.5,
    alpha = 0.75
  ) +
  ylab("Mean CRPS") +
  theme(aspect.ratio = 1) +
  xlab("estimate SD/error SD")

# .gg_save_as(p2a,size=std_size$third,formats = "pdf")

# ggplot(analysis_data, aes(x=unbiased_percent_95_ci, y=unbiased_mean_crps, colour=abs(bias)))+geom_point()+.gg_scale_x_logit()+
#   facet_wrap(~sd_label)+scale_y_log10()
#
# ggplot(analysis_data, aes(x=directed_was, y=mean_crps, colour=abs(bias)))+geom_point()+
#   facet_wrap(~sd_label)+scale_y_log10()
#
# ggplot(analysis_data, aes(x=conf_vs_error, y=mean_crps, colour=abs(bias)))+geom_point()+
#   facet_wrap(~sd_label)+scale_y_log10()+scale_x_log10()
#
# ggplot(analysis_data, aes(x=conf_vs_error, y=unbiased_mean_crps, colour=confidence))+geom_point()+scale_x_log10()+
#   facet_wrap(~sd_label)+scale_y_log10()
#
# ggplot(analysis_data, aes(x=directed_was, y=unbiased_mean_crps, colour=confidence))+geom_point()+
#   facet_wrap(~sd_label)+scale_y_log10()

p3b = ggplot(
  analysis_data,
  aes(x = conf_vs_error, y = directed_was, colour = conf_sd_cv)
) +
  geom_point() +
  scale_x_log10() +
  # facet_wrap(~sd_label)+
  scale_color_viridis_c(name = "CI variation") +
  xlab("estimate SD/error SD") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1) +
  annotate(
    geom = "label",
    label = "overconfident >",
    x = 1,
    y = 0,
    hjust = -0.5,
    vjust = 0.5,
    angle = 90,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< conservative",
    x = 1,
    y = 0,
    hjust = 1.5,
    vjust = 0.5,
    angle = 90,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = 0,
    hjust = -0.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = 0,
    hjust = 1.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  theme(aspect.ratio = 1) +
  ylab("directed PIT wasserstein")

p3a = ggplot(
  analysis_data,
  aes(x = conf_vs_error, y = unbiased_was, colour = conf_sd_cv)
) +
  geom_point() +
  scale_x_log10() +
  # facet_wrap(~sd_label)+
  scale_color_viridis_c(name = "CI variation") +
  xlab("estimate SD/error SD") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 1) +
  annotate(
    geom = "label",
    label = "< better",
    x = 1,
    y = 0,
    hjust = -1.5,
    vjust = 0.5,
    angle = 90,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = 0,
    hjust = -0.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = 0,
    hjust = 1.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  theme(aspect.ratio = 1) +
  ylab("bias adjusted\nPIT wasserstein")

# ggplot(analysis_data, aes(x=unbiased_percent_95_ci, y=directed_was, colour=conf_sd_cv))+
#   geom_point()+
#   facet_wrap(~sd_label)+
#   geom_vline(xintercept = 0.95)+geom_hline(yintercept=0)+
#   annotate(geom="label",label="overconfident \u2192", x=0.95,y=0,hjust=-0.6,vjust=0.5,angle=90,alpha=0.75)+
#   annotate(geom="label",label="\u2190 conservative", x=0.95,y=0,hjust=1.6,vjust=0.5,angle=90,alpha=0.75)+
#   annotate(geom="label",label="conservative \u2192", x=0.95,y=0,hjust=-0.6,vjust=0.5,alpha=0.75)+
#   annotate(geom="label",label="\u2190 overconfident", x=0.95,y=0,hjust=1.6,vjust=0.5,alpha=0.75)+
#   xlab("bias adjusted\ncoverage probability")+
#   .gg_scale_x_logit() + theme(aspect.ratio=1)+
#   scale_color_viridis_c(name="CI variation")+
#   ylab("directed PIT wasserstein")

p2b = ggplot(
  analysis_data,
  aes(y = percent_50_ci, x = conf_vs_error, colour = conf_sd_cv)
) +
  geom_point() +
  # facet_wrap(~sd_label)+
  geom_vline(xintercept = 1) +
  geom_hline(yintercept = 0.5) +
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = 0.5,
    hjust = -0.5,
    vjust = 0.5,
    angle = 90,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = 0.5,
    hjust = 1.5,
    vjust = 0.5,
    angle = 90,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = 0.5,
    hjust = -0.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = 0.5,
    hjust = 1.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  ylab("50% coverage probability") +
  # .gg_scale_y_logit(sf=3) +
  theme(aspect.ratio = 1) +
  scale_color_viridis_c(name = "CI variation") +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("estimate SD/error SD")

# p2 = p2a+p2b+patchwork::plot_layout(ncol=1, guides = "collect",axes="collect")+patchwork::plot_annotation(tag_levels = "A")

p2 = p2a +
  p2b +
  p3a +
  p3b +
  patchwork::plot_layout(ncol = 2, guides = "collect", axes = "collect") +
  patchwork::plot_annotation(tag_levels = "A")

.gg_save_as(
  p2,
  here::here("s2/fig/fig3-metric-compare"),
  size = std_size$two_third,
  formats = c("pdf", "eps")
)


# p3 = p3a +
#   p3b +
#   patchwork::plot_layout(ncol = 1, guides = "collect", axes = "collect") +
#   patchwork::plot_annotation(tag_levels = "A")
# .gg_save_as(
#   p3,
#   here::here("s2/fig/fig4-metric-compare"),
#   size = std_size$two_third,
#   formats = c("pdf", "eps")
# )

# $$
#   \epsilon_i = \overline{x_i} - \int_0^\infty{x f_i(x) dx}
# $$

# ggplot(analysis_data, aes(x=ifelse(unbiased_percent_95_ci==1,1-1/length(true_value),unbiased_percent_95_ci), y=directed_was, colour=biased))+
#   geom_point()+
#   facet_wrap(~sd_label)+
#   geom_vline(xintercept = 0.95)+geom_hline(yintercept=0)+
#   annotate(geom="label",label="overconfident \u2192", x=0.95,y=0,hjust=-0.6,vjust=0.5,angle=90,alpha=0.75)+
#   annotate(geom="label",label="\u2190 convervative", x=0.95,y=0,hjust=1.6,vjust=0.5,angle=90,alpha=0.75)+
#   annotate(geom="label",label="convervative \u2192", x=0.95,y=0,hjust=-0.6,vjust=0.5,alpha=0.75)+
#   annotate(geom="label",label="\u2190 overconfident", x=0.95,y=0,hjust=1.6,vjust=0.5,alpha=0.75)+
#   xlab("bias adjusted\ncoverage probability")+
#   .gg_scale_x_logit() + theme(aspect.ratio=1)+
#   ylab("directed PIT wasserstein")
#
# ggplot(analysis_data, aes(x=ifelse(unbiased_percent_95_ci==1,1-1/length(true_value),unbiased_percent_95_ci), y=directed_was, colour=confidence))+
#   geom_point()+
#   facet_wrap(~sd_label)+
#   geom_vline(xintercept = 0.95)+geom_hline(yintercept=0)+
#   annotate(geom="label",label="overconfident \u2192", x=0.95,y=0,hjust=-0.6,vjust=0.5,angle=90,alpha=0.75)+
#   annotate(geom="label",label="\u2190 convervative", x=0.95,y=0,hjust=1.6,vjust=0.5,angle=90,alpha=0.75)+
#   annotate(geom="label",label="convervative \u2192", x=0.95,y=0,hjust=-0.6,vjust=0.5,alpha=0.75)+
#   annotate(geom="label",label="\u2190 overconfident", x=0.95,y=0,hjust=1.6,vjust=0.5,alpha=0.75)+
#   xlab("bias adjusted\ncoverage probability")+
#   .gg_scale_x_logit() + theme(aspect.ratio=1)+
#   ylab("directed PIT wasserstein")

# relationship between
# analysis_data %>%
#   group_by(sd_label) %>%
#   transmute(
#     conf_factor = log(est_error_sd / conf_sd_median),
#     directed_was,
#     wass_factor = conf_factor / directed_was
#   ) %>%
#   filter(conf_factor > 0) %>%
#   summarise(
#     wass_factor_mean = mean(wass_factor),
#     wass_factor_sd = sd(wass_factor)
#   )

# Figure 4 ----

no_bias_scenarios = tibble(
  sd_sd_mean = c(0),
  sd_label = forcats::as_factor(c("uniform CIs"))
) %>%
  group_by(sd_label, sd_sd_mean) %>%
  reframe(
    est_error_mean = 0,
    est_error_sd = exp(runif(1000, -log(conf_range), log(conf_range))),
    conf_sd_median = exp(runif(1000, -log(conf_range), log(conf_range))),
    conf_sd_cv = runif(1000, 0, sd_sd_mean)
  ) %>%
  ungroup()
no_bias_data = run_scenario(true_value, no_bias_scenarios) %>%
  mutate(
    # convert mean to median assuming uninformed beta distribution prior to
    # allow values with zero counts to be represented
    unbiased_percent_50_ci = (unbiased_percent_50_ci * n + 1 / 3) / (n + 2 / 3),
    conf_vs_error = conf_sd_median / est_error_sd
  )

p4b = ggplot(
  no_bias_data,
  aes(x = conf_vs_error, y = mean_crps, colour = est_error_sd)
) +
  geom_point() +
  scale_x_log10() +
  # facet_wrap(~sd_label)+
  geom_vline(xintercept = 1) +
  scale_color_viridis_c(option = "magma", name = "Error SD") +
  # ggpp::annotate(geom=ggpp::GeomLabelNpc, label="conservative \u2192", npcx=0.6,npcy=0.95,hjust="left",vjust="center",alpha=0.75, size=6, size.unit="pt", label.r=unit(0,"pt"))+
  # ggpp::annotate(geom=ggpp::GeomLabelNpc, label="\u2190 overconfident", npcx=0.4,npcy=0.95,hjust="right",vjust="center",alpha=0.75, size=6, size.unit="pt", label.r=unit(0,"pt"))+
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = Inf,
    hjust = -0.5,
    vjust = 1.5,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = Inf,
    hjust = 1.5,
    vjust = 1.5,
    alpha = 0.75
  ) +
  ylab("Mean CRPS") +
  theme(aspect.ratio = 1) +
  xlab("estimate SD/error SD")

p4a = ggplot(
  no_bias_data,
  aes(y = percent_50_ci, x = conf_vs_error, colour = est_error_sd)
) +
  geom_point() +
  # facet_wrap(~sd_label)+
  geom_vline(xintercept = 1) +
  geom_hline(yintercept = 0.5) +
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = 0.5,
    hjust = -0.5,
    vjust = 0.5,
    angle = 90,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = 0.5,
    hjust = 1.5,
    vjust = 0.5,
    angle = 90,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "conservative >",
    x = 1,
    y = 0.5,
    hjust = -0.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  annotate(
    geom = "label",
    label = "< overconfident",
    x = 1,
    y = 0.5,
    hjust = 1.5,
    vjust = 0.5,
    alpha = 0.75
  ) +
  ylab("50% coverage probability") +
  scale_color_viridis_c(option = "magma", name = "Error SD") +
  # .gg_scale_y_logit(sf=3) +
  theme(aspect.ratio = 1) +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("estimate SD/error SD")

p4 = p4b +
  p4a +
  patchwork::plot_layout(ncol = 2, guides = "collect", axes = "collect") +
  patchwork::plot_annotation(tag_levels = "A")

.gg_save_as(
  p4,
  here::here("s2/fig/fig4-metric-no-bias"),
  size = std_size$third,
  formats = c("pdf", "eps")
)
# Figure 1 ----

pdata = tibble(
  x = seq(-1, 5, 0.01),
  y = dnorm(x, 2, 0.5)
) %>%
  mutate(y = y / max(y))

lndata = tibble(
  x = seq(0, 2, 0.01),
  y = dlnorm(x, log(0.35), 0.25)
) %>%
  mutate(y = y / max(y))

pest = tibble(
  x = seq(-1, 5, 0.01),
  y = dnorm(x, 2.5, 0.25)
) %>%
  mutate(y = y / max(y))

pBack = ggplot() +
  coord_cartesian(ylim = c(-0.25, NA), xlim = c(0.5, 4.5)) +
  theme_void() +
  # True value
  geom_vline(xintercept = 1, colour = "blue", size = 0.75) +
  annotate(
    geom = "text",
    x = 1,
    y = 0.075,
    label = "true value",
    hjust = 0,
    vjust = -0.5,
    angle = 90,
    colour = "blue"
  ) +

  # Spread distribution
  geom_ribbon(
    data = lndata,
    aes(x = x + 2.5, ymin = -0.1, ymax = -0.1 - y * 0.1),
    fill = "orange",
    colour = "orange",
    alpha = 0.2
  ) +
  annotate(
    geom = "text",
    x = 3,
    y = -0.15,
    label = "estimator\nspread distribution",
    hjust = 0,
    vjust = 0.5,
    colour = "orange"
  ) +
  # Error Function
  geom_ribbon(
    data = pdata,
    aes(x = x, ymax = y * 0.5),
    ymin = 0,
    fill = "red",
    colour = "red",
    alpha = 0.2
  ) +
  geom_segment(
    data = pdata %>% filter(x == 2),
    aes(x = x, y = y * 0.5, xend = x),
    yend = 0,
    colour = "red",
    linetype = "dashed"
  ) +
  # annotate(geom = "segment", x=1,xend=3,y=0,yend=0, arrow=arrow(angle=90,length=unit(2,"mm"),type="open",ends = "both"),size=0.75,colour="red")+
  annotate(geom = "point", x = 2, y = 0, size = 2, colour = "red") +
  annotate(
    geom = "text",
    x = 2.5,
    y = 0.75 * 0.5,
    label = "estimator\nerror function",
    hjust = 0,
    vjust = 0.5,
    colour = "red"
  ) +

  # Estimate example
  geom_ribbon(
    data = pest,
    aes(x = x, ymin = -0.05, ymax = -0.05 + y * 0.2),
    fill = "black",
    colour = "black",
    alpha = 0.2
  ) +
  annotate(
    geom = "segment",
    x = 2,
    xend = 3,
    y = -0.05,
    yend = -0.05,
    arrow = arrow(
      angle = 90,
      length = unit(2, "mm"),
      type = "open",
      ends = "both"
    ),
    size = 0.75,
    colour = "black"
  ) +
  annotate(geom = "point", x = 2.5, y = -0.05, size = 2) +
  annotate(
    geom = "text",
    x = 3.75,
    y = -0.05,
    label = "estimate + CIs",
    hjust = 0,
    vjust = 1.5
  ) +

  # Error linking arrow
  geom_segment(
    data = pdata %>% filter(x == 2.5),
    aes(x = x, y = y * 0.5, xend = x),
    yend = -0.05,
    arrow = arrow(angle = 20, length = unit(2, "mm"), type = "closed"),
    colour = "grey20"
  ) +
  geom_point(
    data = pdata %>% filter(x == 2.5),
    aes(x = x, y = y * 0.5),
    size = 1
  ) +

  # CI linking arrow
  geom_segment(
    data = lndata %>% filter(x == 2.75 - 2.5),
    aes(x = x + 2.5, xend = x + 2.5, y = -0.1 - y * 0.1),
    yend = -0.065,
    arrow = arrow(angle = 20, length = unit(2, "mm"), type = "closed"),
    colour = "grey20"
  ) +
  geom_point(
    data = lndata %>% filter(x == 2.75 - 2.5),
    aes(x = x + 2.5, y = -0.1 - y * 0.1),
    size = 1
  ) +
  annotate(
    geom = "segment",
    x = 2.5,
    xend = 2.75,
    y = -0.065,
    yend = -0.065,
    arrow = arrow(
      angle = 20,
      length = unit(2, "mm"),
      type = "closed",
      ends = "both"
    ),
    colour = "grey20"
  ) +

  # Bias arrow
  annotate(
    geom = "segment",
    x = 1,
    xend = 2,
    y = 0.25,
    yend = 0.25,
    arrow = arrow(
      angle = 20,
      length = unit(2, "mm"),
      type = "closed",
      ends = "both"
    ),
    colour = "grey20"
  ) +
  annotate(
    geom = "text",
    x = 1.5,
    y = 0.25,
    label = "bias",
    hjust = 0.5,
    vjust = -1.25
  )
pBack
.gg_save_as(
  pBack,
  here::here("s2/fig/fig1-introduction"),
  size = std_size$sixth,
  formats = c("pdf", "eps")
)

# A data generating function G(x)
# A estimator F(x)

# Assume G(x) = F(x-bias) - i.e. F is biased but not under or over dispersed
# implies G'(y) = F'(y) - bias

# PIT is F(y) where y is a sample from G
# PIT histogram distribution:
# F(G'(y)-bias) where y \sim Unif(0,1)
# F(F'(y) - bias)
