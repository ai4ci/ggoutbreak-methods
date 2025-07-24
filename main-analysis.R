library(tidyverse)
library(patchwork)
# devtools::load_all("~/Git/interfacer")
#devtools::load_all("~/Git/ggoutbreak")
library(ggoutbreak)

future::plan(future::multisession, workers = 12)

source(here::here("R/utils.R"))

.gg_pedantic(
  size = 0.5,
  lineSize = 0.25,
  fontSize = 6,
  panel.grid = element_blank()
)

# Infectivity profile ----
# generate 100 epiestim compatible generation intervals with a mean of (approx) 5 and a SD of
# (approx) 3. The average of these will be the `truth`
withr::with_seed(200, {
  ip = ggoutbreak::make_gamma_ip(
    median_of_mean = 5,
    lower_ci_of_mean = 4,
    upper_ci_of_mean = 6,
    median_of_sd = 2,
    lower_ci_of_sd = 1,
    upper_ci_of_sd = 3,
    n_boots = 100,
    epiestim_compat = TRUE
  )
})

sfig1 = ggoutbreak::plot_ip(ip, alpha = 0.2) +
  xlab("days since primary infection (τ)") +
  scale_x_continuous(breaks = .gg_integer_breaks())


.gg_save_as(
  sfig1,
  filename = here::here("s1/fig/fig1-simulated-inf-prof"),
  size = std_size$third,
  formats = c("pdf", "eps")
)

# Reproduction number as parameters ----

# Figure 1: Qualitative comparison of estimation methods ----
# An example time series with ascertainment issues plot of cases as row 1.

# generate a demonstration R_t time-series as input parameters
changes = tibble::tibble(
  t = ggoutbreak::as.time_period(c(0, 20, 40, 60, 80, 110), "1 day"),
  R = c(1.5, 1.3, 0.9, 1.2, 0.8, 1.2)
)

# A branching process model with the
bpm = sim_branching_process(
  changes = changes,
  max_time = 120,
  fn_ip = ~ip,
  seed = 100,
  fn_imports = ~ ifelse(.x == 0, 30, 0)
)
events = attr(bpm, "events")
count_bpm = bpm %>% sim_summarise_linelist()

kappas = c(low = 0.1, medium = 0.4, high = 0.7)

# # Dispersion parameters coefficient of variation
# This is just for reference really so that we can report something easily.
# lapply(kappas, \(kappa) {
#   tmp = ggoutbreak:::.reparam_beta(0.7,kappa)
#   var = (tmp$alpha*tmp$beta)/((tmp$alpha+tmp$beta)^2*(tmp$alpha+tmp$beta+1))
#   mean = tmp$alpha/(tmp$alpha+tmp$beta)
#   format( sqrt(var)/(mean), digits = 4,scientific = TRUE)
# })

withr::with_seed(101, {
  obs_bpm = bind_rows(lapply(seq_along(kappas), function(i) {
    # generate observed data for each of the ascertainment values
    # This is done by applying a per day random factor from a beta
    # distribution with mean 0.7 and kappa value as defined by `kappas`
    kappa = kappas[i]
    count_bpm %>%
      sim_apply_ascertainment(~ rbeta2(.x, 0.7, kappa)) %>%
      mutate(variation = names(kappas)[i])
  })) %>%
    group_by(variation) %>%
    mutate(
      variation = factor(variation, levels = names(kappas))
    )

  # Generate an estimate of incidence from the observed BPM with noise
  # using a poisson model with log link function
  model_bpm = obs_bpm %>%
    poisson_locfit_model(window = 14, deg = 2) %>%
    mutate(model = "Incidence model")

  # calculate the RT from incidence
  rt_incidence_bpm = model_bpm %>%
    rt_from_incidence(ip = ip, approx = FALSE) %>%
    mutate(model = "Rₜ from incidence")

  # calculate the RT using epiestim
  epiestim_bpm = obs_bpm %>%
    rt_epiestim(ip = ip, window = 14) %>%
    mutate(model = "EpiEstim")
  rt_comparison = bind_rows(rt_incidence_bpm, epiestim_bpm) %>%
    group_by(model, variation)
})

p1 = plot_incidence(model_bpm, raw = obs_bpm, events = events) +
  facet_grid(model ~ variation)
p2 = plot_rt(rt_comparison, events = events) +
  facet_grid(model ~ variation) +
  geom_step(
    data = obs_bpm,
    mapping = aes(x = as.Date(time), y = rt.weighted),
    colour = "red"
  ) +
  coord_cartesian(ylim = c(0.6, 2))

fig1 =
  p1 +
  theme(axis.text.x.bottom = element_blank()) +
  p2 +
  theme(axis.text.x.top = element_blank()) +
  patchwork::plot_layout(ncol = 1, axes = "collect", heights = c(1, 2)) +
  patchwork::plot_annotation(tag_levels = "A")

.gg_save_as(
  fig1,
  filename = here::here("main/fig/fig1-noise-qualitative"),
  size = std_size$half,
  formats = c("pdf", "eps")
)

# ylab("simulation Rₜ")+xlab(NULL)+
# scale_x_date(date_breaks = "2 week",date_labels = "%d %b")+

# Figure 2: Quantitative comparisons ----

## generate 5 outbreak scenarios ----

# Setup caches and cached functions
# cache = memoise::cache_filesystem(here::here("cache"))

setup_scenarios = function(n, seed) {
  withr::with_seed(seed, {
    tibble(scenario = 1:n) %>%
      mutate(
        changes = purrr::map(
          scenario,
          ~ tibble::tibble(
            # gap between changes in R_t, 4 entries each scenario
            t = cumsum(c(0, rpois(3, 20))),
            # Rt values at change. oscillates between positive and negative
            R = exp(runif(4, 0, log(1.5)) * ((-1)^(0:3)))
          )
        )
      ) %>%
      tidyr::crossing(seed = runif(50))
  })
}

#setup_scen = memoise::memoise(setup_scenarios, cache = cache)
#sim_bpm = memoise::memoise(sim_branching_process, cache = cache)

setup_scen = setup_scenarios
sim_bpm = sim_branching_process


qual_df = setup_scen(5, 100)
max_time = max(c(80, purrr::map_dbl(qual_df$changes, ~ max(.x$t))))

## build the simulations ---
## This is 5 random scenarios with 50 random seeds (250 replicas)
## qual_df2 is a nested DF with 250 rows once per replica
## nested `bpm` value is a individual level line list
qual_df2 = qual_df %>%
  dplyr::mutate(
    bpm = furrr::future_map2(
      changes,
      seed,
      ~ suppressMessages(sim_bpm(
        changes = .x,
        seed = .y,
        max_time = max_time,
        fn_ip = ~ip,
        fn_imports = ~ ifelse(.x == 0, 30, 0)
      )),
      .progress = TRUE
    )
  )

# Summarise the BPM linelists to a daily count
qual_df2 = qual_df2 %>%
  dplyr::mutate(
    summ_bpm = purrr::map(bpm, ~ sim_summarise_linelist(.x), .progress = TRUE)
  )

options("ggoutbreak.keep_cdf" = FALSE)

## generate a noisy observed set ----
# The 250 bpm daily counts have 2 types of random noise applied to them
# based on kappas as before. This will be 750.

qual_df3 = qual_df2 %>%
  dplyr::select(-bpm) %>%
  tidyr::crossing(asc_kappa = c(0.1, 0.4, 0.7))

qual_df4 = qual_df3 %>%
  dplyr::mutate(
    asc = purrr::pmap(., \(summ_bpm, asc_kappa, seed, ...) {
      withr::with_seed(seed, rbeta2(nrow(summ_bpm), 0.7, asc_kappa))
    }),
    obs_bpm = purrr::map2(
      summ_bpm,
      asc,
      ~ sim_apply_ascertainment(.x, \(t) .y),
      .progress = TRUE
    )
  )

## estimate and compare scores for different methods for all 750 scenarios ----

do_comparison = function(scenarios, window, use_lags, date_cutoff = 20) {
  # window = 14
  # scenarios = qual_df4 %>% filter(scenario==1, seed==min(seed))

  if (use_lags) {
    incid_lag = quantify_lag(
      ~ poisson_locfit_model(.x, window = window, deg = 2) %>%
        rt_from_incidence(ip = .y, approx = FALSE),
      ip = ip
    )
    # We are using ggoutbreak reimplementation of CORI method here as is faster.
    # This is shown to be identical to EpiEstim when not using random resampling
    # to combine estimates.
    cori_lag = quantify_lag(
      ~ rt_cori(.x, ip = .y, window = window, epiestim_compat = TRUE),
      ip = ip
    )
  } else {
    incid_lag = NULL
    cori_lag = NULL
  }

  # qual_df4 %>% glimpse()
  qual_df5 = scenarios %>%
    dplyr::mutate(
      model_bpm = purrr::map(
        obs_bpm,
        ~ poisson_locfit_model(.x, window = window, deg = 2),
        .progress = "Modelling incidence 1/3"
      )
    )
  message("Rt from incidence (+scoring) 2/3")
  # Incidence model
  qual_df6 = qual_df5 %>%
    dplyr::mutate(
      rt_incidence_bpm = furrr::future_map2(
        model_bpm,
        summ_bpm,
        ~ {
          withr::with_options(list("ggoutbreak.keep_cdf" = TRUE), {
            tmp = ggoutbreak::rt_from_incidence(.x, ip = ip, approx = TRUE)
            tmp2 = ggoutbreak::score_estimate(
              est = tmp |> dplyr::filter(time > date_cutoff),
              obs = .y |> dplyr::rename(rt.obs = rt.weighted),
              lags = incid_lag,
              # we have 50 replicates so we don't need lots of boots
              # we are not summarising for this estimate as we will combine it
              # with other replicates, in the plot
              bootstraps = 20,
              raw_bootstraps = TRUE
            )
            tmp = tmp |> dplyr::select(-tidyselect::ends_with(".cdf"))
          })
          return(list(est = tmp, score = tmp2))
        },
        .progress = TRUE #"Rt from incidence (+scoring) 2/3"
      )
    )

  # Cori method model - equivalent to EpiEstim with Quantile Bias score.
  message("EpiEstim (+scoring) 3/3")
  qual_df7 = qual_df6 %>%
    dplyr::mutate(
      rt_cori_bpm = furrr::future_map2(
        obs_bpm,
        summ_bpm,
        ~ {
          withr::with_options(list("ggoutbreak.keep_cdf" = TRUE), {
            tmp = ggoutbreak::rt_cori(
              .x,
              ip = ip,
              window = window,
              epiestim_compat = TRUE
            )
            tmp2 = ggoutbreak::score_estimate(
              est = tmp %>% filter(time > date_cutoff),
              obs = .y %>% rename(rt.obs = rt.weighted),
              lags = cori_lag,
              # we have 50 replicates so we don't need lots of boots
              # we are not summarising for this estimate as we will combine it
              # with other replicates, in the plot
              bootstraps = 20,
              raw_bootstraps = TRUE
            )
            tmp = tmp %>% select(-tidyselect::ends_with(".cdf"))
          })
          return(list(est = tmp, score = tmp2))
        },
        .progress = TRUE # "EpiEstim (+scoring) 3/3"
      )
    )

  tmp1 = qual_df7 %>%
    mutate(
      rt_bpm = purrr::map(rt_incidence_bpm, ~ .x$est),
      rt_score = purrr::map(
        rt_incidence_bpm,
        ~ .x$score %>% filter(.type == "rt")
      ),
      method = "Rₜ from incidence"
    ) %>%
    select(
      -rt_cori_bpm,
      -rt_incidence_bpm #, -rt_epi_bpm
    )

  tmp2 = qual_df7 %>%
    mutate(
      rt_bpm = purrr::map(rt_cori_bpm, ~ .x$est),
      rt_score = purrr::map(
        rt_cori_bpm,
        ~ .x$score %>% filter(.type == "rt")
      ),
      method = "EpiEstim"
    ) %>%
    select(
      -rt_cori_bpm,
      -rt_incidence_bpm #, -rt_epi_bpm
    ) %>%
    mutate(
      seed_id = min_rank(seed)
    )

  # qual_df9 = dplyr::bind_rows(tmp1, tmp2)

  return(dplyr::bind_rows(tmp1, tmp2))
}


if (interactive()) {
  qual_df9 = do_comparison(qual_df4, 14, TRUE)
  saveRDS(qual_df9, here::here("cache/full_estimates_and_scoring.Rds"))
} else {
  qual_df9 = readRDS(here::here("cache/full_estimates_and_scoring.Rds"))
}

if (interactive()) {
  qual_df9_7_day = do_comparison(qual_df4, 7, TRUE)
  saveRDS(qual_df9_7_day, here::here("cache/7_day_estimates_and_scoring.Rds"))
} else {
  qual_df9_7_day = readRDS(here::here("cache/7_day_estimates_and_scoring.Rds"))
}

if (interactive()) {
  qual_df9_no_lag = do_comparison(qual_df4, 14, FALSE)
  saveRDS(qual_df9_no_lag, here::here("cache/no_lag_estimates_and_scoring.Rds"))
} else {
  qual_df9_no_lag = readRDS(here::here(
    "cache/no_lag_estimates_and_scoring.Rds"
  ))
}


incid_lag = quantify_lag(
  ~ poisson_locfit_model(.x, window = 14, deg = 2) %>%
    rt_from_incidence(ip = .y, approx = FALSE),
  ip = ip
)
cori_lag = quantify_lag(
  ~ rt_cori(.x, ip = .y, window = 14, epiestim_compat = TRUE),
  ip = ip
)
incid_lag_2 = quantify_lag(
  ~ poisson_locfit_model(.x, window = 7, deg = 2) %>%
    rt_from_incidence(ip = .y, approx = FALSE),
  ip = ip
)
cori_lag_2 = quantify_lag(
  ~ rt_cori(.x, ip = .y, window = 7, epiestim_compat = TRUE),
  ip = ip
)

lagtable = bind_rows(
  incid_lag %>% mutate(method = "Rₜ from incidence", window = 14),
  cori_lag %>% mutate(method = "EpiEstim", window = 14),
  incid_lag_2 %>% mutate(method = "Rₜ from incidence", window = 7),
  cori_lag_2 %>% mutate(method = "EpiEstim", window = 7)
) %>%
  filter(estimate == "rt") %>%
  ungroup() %>%
  select(-estimate)

lagtable %>%
  mutate(window = sprintf("%d", window), label = "Window (days)") %>%
  .hux_tidy(c("method"), c("label", "window")) %>%
  .hux_save_as("main/fig/tab1-lags.pdf", size = std_size$quarter_portrait)

lagplot = bind_rows(
  attributes(incid_lag)$data %>%
    mutate(method = "Rₜ from incidence", window = 14),
  attributes(cori_lag)$data %>%
    mutate(method = "EpiEstim", window = 14),
  attributes(incid_lag_2)$data %>%
    mutate(method = "Rₜ from incidence", window = 7),
  attributes(cori_lag_2)$data %>%
    mutate(method = "EpiEstim", window = 7)
) %>%
  filter(.type == "rt") %>%
  ggplot(aes(x = lag, y = lagged_rmse, colour = as.factor(window))) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_color_brewer(palette = "Dark2", name = "window") +
  facet_wrap(~method) +
  ylab("Root mean square error") +
  xlab("Estimate delay (days)")


.gg_save_as(
  lagplot,
  filename = here::here("s1/fig/fig2-lag-plot"),
  size = std_size$third,
  formats = c("pdf", "eps")
)


# Re-jig nested list format and make long
# tmp3 = qual_df8 %>% mutate(
#   rt_bpm = purrr::map(rt_cori_bpm, ~ .x$est),
#   rt_score = purrr::map(rt_cori_bpm, ~ .x$score %>% filter(estimate=="rt")),
#   method = "EpiEstim"
# ) %>% select(
#   -rt_cori_bpm, -rt_epi_bpm, -rt_incidence_bpm
# )
#
# qual_df9 = dplyr::bind_rows(tmp1,tmp2,tmp3)
#
# qual_df9 = qual_df9 %>% mutate(seed_id = min_rank(seed)) %>% glimpse()
# n_distinct(qual_df9$seed)

.kappa_lbl = function(asc_kappa) {
  matched = sapply(asc_kappa, \(.x) which(kappas == .x))
  return(factor(names(kappas)[matched], levels = c("low", "medium", "high")))
}

do_figure_2 = function(qual_df9) {
  tmp = qual_df9 %>%
    mutate(seed_id = factor(seed, labels = as.character(1:50))) %>%
    select(method, asc_kappa, scenario, seed_id, rt_score) %>%
    unnest(rt_score)

  tmp2 = qual_df9 %>%
    mutate(seed_id = factor(seed, labels = as.character(1:50))) %>%
    select(method, asc_kappa, scenario, seed_id, summ_bpm) %>%
    unnest(summ_bpm)

  tmp_rt = qual_df9 %>%
    mutate(seed_id = factor(seed, labels = as.character(1:50))) %>%
    select(method, asc_kappa, scenario, seed_id, rt_bpm) %>%
    unnest(rt_bpm)

  min_max = tmp_rt %>%
    group_by(scenario, seed_id) %>%
    summarise(rt.mean = mean(rt.0.5)) %>%
    arrange(rt.mean) %>%
    filter(row_number() == 1 | row_number() == n())

  s0 = ggplot(
    tmp_rt %>% semi_join(min_max, by = c("scenario", "seed_id")),
    aes(
      x = time,
      y = rt.0.5,
      group = interaction(asc_kappa, scenario, seed_id, method),
      colour = method
    )
  ) +
    geom_ribbon(
      aes(ymin = rt.0.025, ymax = rt.0.975, fill = method),
      colour = NA,
      alpha = 0.1
    ) +
    geom_line() +
    geom_line(data = tmp2, aes(y = rt.weighted), colour = "black") +
    geom_hline(yintercept = 1, colour = "grey30") +
    facet_grid(scenario ~ .kappa_lbl(asc_kappa)) +
    coord_cartesian(ylim = c(0, 2)) +
    ylab("Rₜ")

  # CRPS for all individual estimates and scenarios combined
  p1 = ggplot(
    tmp,
    aes(y = mean_crps, fill = method, x = .kappa_lbl(asc_kappa))
  ) +
    geom_boxplot(outliers = FALSE) +
    #theme(legend.position = "bottom")+
    xlab("Dispersion") +
    ylab("CRPS") +
    coord_cartesian(ylim = c(0, NA))

  p2 = ggplot(
    tmp,
    aes(
      y = (mean_bias - 1) * 100,
      fill = method,
      x = .kappa_lbl(asc_kappa)
    )
  ) +
    geom_boxplot(outliers = FALSE) +
    #theme(legend.position = "bottom")+
    xlab("Dispersion") +
    geom_hline(yintercept = 0, colour = "grey30") +
    ylab("Proportional bias (%)")

  p3 = ggplot(
    tmp,
    aes(
      y = mean_prediction_interval_width_50,
      fill = method,
      x = .kappa_lbl(asc_kappa)
    )
  ) +
    geom_boxplot(outliers = FALSE) +
    #theme(legend.position = "bottom")+
    xlab("Dispersion") +
    coord_cartesian(ylim = c(0, NA)) +
    ylab("Prediction interval width")

  p4 = ggplot(
    tmp,
    aes(
      y = unbiased_percent_iqr_coverage,
      fill = method,
      x = .kappa_lbl(asc_kappa)
    )
  ) +
    geom_boxplot(outliers = FALSE) +
    #theme(legend.position = "bottom")+
    xlab("Dispersion") +
    geom_hline(yintercept = 0.5, colour = "grey30") +
    ylab("50% coverage (adj)")

  # p4 = ggplot(tmp, aes(y = directed_pit_was, fill = method, x = .kappa_lbl(asc_kappa))) +
  #   geom_boxplot(outliers = FALSE) +
  #   #theme(legend.position = "bottom")+
  #   xlab("Dispersion") +
  #   geom_hline(yintercept = 0,colour="grey30") +
  #   ylab("Directed PIT Wasserstein (adj)")

  p5 = ggplot(
    tmp,
    aes(y = unbiased_pit_was, fill = method, x = .kappa_lbl(asc_kappa))
  ) +
    geom_boxplot(outliers = FALSE) +
    #theme(legend.position = "bottom")+
    xlab("Dispersion") +
    ylab("PIT Wasserstein (adj)") +
    coord_cartesian(ylim = c(0, NA))

  p6 = ggplot(
    tmp,
    aes(
      y = threshold_misclassification_probability,
      fill = method,
      x = .kappa_lbl(asc_kappa)
    )
  ) +
    geom_boxplot(outliers = FALSE) +
    #theme(legend.position = "bottom")+
    xlab("Dispersion") +
    ylab("Threshold misclassification") +
    coord_cartesian(ylim = c(0, NA))

  pout = p1 +
    guides(fill = guide_none()) +
    p2 +
    # p2a +
    p3 +
    p4 +
    p5 +
    p6 +
    patchwork::plot_layout(guides = "collect", axes = "collect", ncol = 3) +
    plot_annotation(tag_levels = "A")

  # pout

  # aggregate CPRS for scenarios and methods and p-values
  # comp = tmp %>%
  #   group_by(asc_kappa, scenario) %>%
  #   group_modify(function(d, g, ...) {
  #     if (n_distinct(d$method) == 2) {
  #       tmp = t.test(mean_crps ~ method, d)
  #       return(broom::tidy(tmp))
  #     } else {
  #       return(tibble::tibble())
  #     }
  #   })
  #
  # comp2 = tmp %>%
  #   group_by(method, asc_kappa, scenario) %>%
  #   summarise(
  #     crps = mean(crps)
  #   )

  # sp0 = tmp %>%
  #   group_by(method, asc_kappa, scenario, seed_id) %>%
  #   reframe(
  #     bias = sort(bias),
  #     uniform = seq(-1, 1, length.out = n())
  #   ) %>%
  #   ggplot(aes(
  #     colour = method,
  #     group = interaction(scenario, seed_id),
  #     x = uniform,
  #     y = bias
  #   )) +
  #   geom_line() +
  #   facet_grid(method ~ .kappa_lbl(asc_kappa)) +
  #   geom_abline(colour = "grey40")

  sp1 = ggplot(
    tmp,
    aes(x = as.factor(scenario), y = mean_crps, fill = method)
  ) +
    geom_boxplot(
      width = 0.5,
      position = position_dodge(width = 0.6),
      outliers = FALSE
    ) +
    ylab("CRPS") +
    xlab("Scenario") +
    facet_wrap(~ .kappa_lbl(asc_kappa))

  sp2 = ggplot(
    tmp,
    aes(
      x = as.factor(scenario),
      y = (mean_bias - 1) * 100,
      fill = method
    )
  ) +
    geom_boxplot(
      width = 0.5,
      position = position_dodge(width = 0.6),
      outliers = FALSE
    ) +
    ylab("% Bias") +
    geom_hline(yintercept = 0, colour = "grey30") +
    xlab("Scenario") +
    facet_wrap(~ .kappa_lbl(asc_kappa))

  sp3 = ggplot(
    tmp,
    aes(
      x = as.factor(scenario),
      y = mean_prediction_interval_width_50,
      fill = method
    )
  ) +
    geom_boxplot(
      width = 0.5,
      position = position_dodge(width = 0.6),
      outliers = FALSE
    ) +
    coord_cartesian(ylim = c(0, NA)) +
    ylab("Prediction interval width") +
    xlab("Scenario") +
    facet_wrap(~ .kappa_lbl(asc_kappa))

  sp4 = ggplot(
    tmp,
    aes(
      x = as.factor(scenario),
      y = percent_iqr_coverage,
      fill = method
    )
  ) +
    geom_boxplot(
      width = 0.5,
      position = position_dodge(width = 0.6),
      outliers = FALSE
    ) +
    ylab("50% coverage") +
    xlab("Scenario") +
    geom_hline(yintercept = 0.5, colour = "grey30") +
    facet_wrap(~ .kappa_lbl(asc_kappa))

  sp5 = ggplot(
    tmp,
    aes(x = as.factor(scenario), y = unbiased_pit_was, fill = method)
  ) +
    geom_boxplot(
      width = 0.5,
      position = position_dodge(width = 0.6),
      outliers = FALSE
    ) +
    ylab("PIT Wasserstein (adj)") +
    xlab("Scenario") +
    facet_wrap(~ .kappa_lbl(asc_kappa))

  sp6 = ggplot(
    tmp,
    aes(
      x = as.factor(scenario),
      y = threshold_misclassification_probability,
      fill = method
    )
  ) +
    geom_boxplot(
      width = 0.5,
      position = position_dodge(width = 0.6),
      outliers = FALSE
    ) +
    #theme(legend.position = "bottom")+
    xlab("Dispersion") +
    ylab("Threshold misclassification") +
    coord_cartesian(ylim = c(0, NA)) +
    facet_wrap(~ .kappa_lbl(asc_kappa))

  spout = sp1 +
    guides(fill = guide_none()) +
    sp2 +
    sp3 +
    sp4 +
    sp5 +
    sp6 +
    patchwork::plot_layout(guides = "collect", axes = "collect", ncol = 1) +
    plot_annotation(tag_levels = "A")

  # # Earth movers distance to uniform from samples
  # emd = function(x, lim = max(abs(x))) {
  #   uniform = seq(-lim, lim, length.out = length(x))
  #
  #   # change < 0 when over confident on positive side
  #   # and > 0 when over confident on negative side
  #   change = uniform - sort(x)
  #
  #   # Divided by 2 here to move from -1 to 1, to 0 to 1 to conform to PIT
  #   # rather than quantile bias.
  #   tibble(
  #     wasserstein = sum(abs(change)) / length(x) / 2,
  #     # centrality will be positive if over confident
  #     centrality = sum(change * sign(uniform)) / length(x) / 2,
  #     n = length(x)
  #   )
  # }
  #
  # # wasserstein_CI = purrr::map(
  # #   seq(2,5,0.1),
  # #   \(n) sapply(1:10000,\(i) emd(runif(10^n,min=-1,max=1),lim=1)),
  # #   .progress = TRUE
  # # )
  # #
  # # wass_lims = tibble(
  # #   n=10^seq(2,5,0.1),
  # #   cis = wasserstein_CI %>% purrr::map(~ quantile(.x,c(0.75,0.95,0.995,0.9995,0.99995)))
  # # )
  # # plot_wass = wass_lims %>% mutate(id=purrr::map(cis, ~ c(0.75,0.95,0.995,0.9995,0.99995))) %>% unnest(c(id,cis))
  # # ggplot(plot_wass,aes(x=n,y=cis,colour=factor(id)))+geom_line()+scale_x_log10()
  #
  # tmp3 = tmp %>%
  #   filter(method %in% c("EpiEstim", "Rₜ from incidence")) %>%
  #   group_by(method, asc_kappa, seed_id) %>%
  #   group_modify(function(d, g, ...) {
  #     out = emd(d$bias, lim = 1)
  #
  #     # KL distance from uniform using locfit.
  #     # lftmp = locfit::locfit(~ locfit::lp(tmpfit,nn = 0.1))
  #     # baseline = integrate(\(x,...) predict(lftmp,newdata=x),lower = -1,upper = 1)
  #     # p_x = baseline$value/2
  #     # KL = integrate(\(x) p_x * log(p_x/predict(lftmp,newdata=x)),lower = -1,upper = 1)$value
  #     #  + integrate(\(x)  predict(lftmp,newdata=x) * log( predict(lftmp,newdata=x) / p_x),lower = -1,upper = 1)$value
  #
  #     return(out)
  #   })
  #
  # # Quantile bias
  #
  # binwidth = 0.05
  #
  # limits = tmp %>%
  #   group_by(method) %>%
  #   count() %>%
  #   mutate(
  #     min_n = qbinom(0.025, size = n, prob = binwidth) / (n * binwidth),
  #     mid_n = qbinom(0.5, size = n, prob = binwidth) / (n * binwidth),
  #     max_n = qbinom(0.975, size = n, prob = binwidth) / (n * binwidth)
  #   ) %>%
  #   ungroup() %>%
  #   summarise(across(where(is.numeric), mean))
  #
  # sp2 = ggplot(tmp, aes(x = mean_quantile_bias, colour = method)) + #,group=interaction(method,seed_id)))+
  #   geom_step(
  #     stat = "bin",
  #     binwidth = 0.05,
  #     mapping = aes(y = after_stat(density)),
  #     direction = "mid"
  #   ) +
  #   xlab("PIT") +
  #   ylab("Density") +
  #   facet_grid(scenario ~ .kappa_lbl(asc_kappa)) +
  #   theme(legend.position = "bottom") +
  #   geom_hline(yintercept = limits$mid_n)
  # # +
  # # geom_hline(yintercept = limits$min_n, linetype="dashed")+
  # # geom_hline(yintercept = limits$max_n, linetype="dashed")
  # #geom_segment(aes(xend=bias,y=-0.1-ifelse(method=="EpiEstim",0.1,0), yend=-0.2-ifelse(method=="EpiEstim",0.1,0)),alpha=0.1)
  #
  # wass_ci = bind_rows(lapply(1:10000, \(i) {
  #   emd(runif(3500, min = -1, max = 1), lim = 1)
  # }))
  #
  # p2 = ggplot(
  #   tmp3,
  #   aes(x = .kappa_lbl(asc_kappa), y = wasserstein, fill = method)
  # ) +
  #   stat_summary(
  #     fun = ~ quantile(.x, 0.5),
  #     geom = "bar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.5,
  #     colour = "black"
  #   ) +
  #   stat_summary(
  #     fun.min = ~ quantile(.x, 0.025),
  #     fun.max = ~ quantile(.x, 0.975),
  #     geom = "errorbar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.1,
  #     colour = "black"
  #   ) +
  #   # theme(legend.position = "bottom")+
  #   geom_hline(
  #     yintercept = quantile(wass_ci$wasserstein, c(0.95)),
  #     linetype = "dotted"
  #   ) +
  #   xlab("Dispersion") +
  #   geom_hline(yintercept = 0) +
  #   ylab("PIT Wasserstein")
  #
  # p3 = ggplot(
  #   tmp3,
  #   aes(x = .kappa_lbl(asc_kappa), y = centrality, fill = method)
  # ) +
  #   stat_summary(
  #     fun = ~ quantile(.x, 0.5),
  #     geom = "bar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.5,
  #     colour = "black"
  #   ) +
  #   stat_summary(
  #     fun.min = ~ quantile(.x, 0.025),
  #     fun.max = ~ quantile(.x, 0.975),
  #     geom = "errorbar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.1,
  #     colour = "black"
  #   ) +
  #   geom_hline(yintercept = 0) +
  #   geom_hline(
  #     yintercept = quantile(wass_ci$centrality, c(0.025)),
  #     linetype = "dotted"
  #   ) +
  #   geom_hline(
  #     yintercept = quantile(wass_ci$centrality, c(0.975)),
  #     linetype = "dotted"
  #   ) +
  #   xlab("Dispersion") +
  #   ylab("PIT Uncertainty")
  #
  # tmp4 = tmp %>%
  #   filter(method %in% c("EpiEstim", "Rₜ from incidence")) %>%
  #   group_by(method, scenario, asc_kappa, seed_id) %>%
  #   group_modify(function(d, g, ...) {
  #     out = emd(d$bias, lim = 1)
  #     return(out)
  #   })
  #
  # wass_ci_2 = bind_rows(lapply(1:10000, \(i) {
  #   emd(runif(3500 / 5, min = -1, max = 1), lim = 1)
  # }))
  #
  # p2a = ggplot(tmp4, aes(x = scenario, y = wasserstein, fill = method)) +
  #   stat_summary(
  #     fun = ~ quantile(.x, 0.5),
  #     geom = "bar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.5,
  #     colour = "black"
  #   ) +
  #   stat_summary(
  #     fun.min = ~ quantile(.x, 0.025),
  #     fun.max = ~ quantile(.x, 0.975),
  #     geom = "errorbar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.1,
  #     colour = "black"
  #   ) +
  #   # theme(legend.position = "bottom")+
  #   geom_hline(
  #     yintercept = quantile(wass_ci_2$wasserstein, c(0.95)),
  #     linetype = "dotted"
  #   ) +
  #   xlab("Scenario") +
  #   geom_hline(yintercept = 0) +
  #   ylab("PIT Wasserstein") +
  #   facet_wrap(~ .kappa_lbl(asc_kappa))
  #
  # p3a = ggplot(tmp4, aes(x = scenario, y = centrality, fill = method)) +
  #   stat_summary(
  #     fun = ~ quantile(.x, 0.5),
  #     geom = "bar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.5,
  #     colour = "black"
  #   ) +
  #   stat_summary(
  #     fun.min = ~ quantile(.x, 0.025),
  #     fun.max = ~ quantile(.x, 0.975),
  #     geom = "errorbar",
  #     position = position_dodge(width = 0.6),
  #     width = 0.1,
  #     colour = "black"
  #   ) +
  #   geom_hline(yintercept = 0) +
  #   geom_hline(
  #     yintercept = quantile(wass_ci_2$centrality, c(0.025)),
  #     linetype = "dotted"
  #   ) +
  #   geom_hline(
  #     yintercept = quantile(wass_ci_2$centrality, c(0.975)),
  #     linetype = "dotted"
  #   ) +
  #   xlab("Scenario") +
  #   ylab("PIT Uncertainty") +
  #   facet_wrap(~ .kappa_lbl(asc_kappa))

  return(list(
    main = pout,
    s0 = s0,
    s1 = spout
    # s2 = sp2,
    # s3 = sp3,
    # p2a = p2a,
    # p3a = p3a
  ))
}

figs = do_figure_2(qual_df9)

# Supplementary
.gg_save_as(
  figs$s0,
  filename = here::here("s1/fig/fig3-scenario-estimates"),
  size = std_size$half,
  formats = c("pdf", "eps")
)
.gg_save_as(
  figs$s1 &
    ggplot2::theme(
      legend.position = "bottom",
      axis.title.y = element_text(size = 6)
    ),
  filename = here::here("s1/fig/fig4-metrics-by-scenario"),
  size = std_size$two_third,
  formats = c("pdf", "eps")
)
# sfigtmp = figs$s0+ylab("Quantile bias")+xlab("Reference")+theme(legend.position = "bottom")+
#   figs$s2+plot_layout(ncol=1,heights = c(3,5))+plot_annotation(tag_levels = "A")

# .gg_save_as(
#   figs$s2 +
#     plot_layout(ncol = 1, heights = c(3, 1)) +
#     plot_annotation(tag_levels = "A"),
#   filename = here::here("s1/fig/fig5-pit-histograms"),
#   size = std_size$two_third,
#   formats = c("pdf", "eps")
# )
#
.gg_save_as(
  figs$main &
    ggplot2::theme(
      legend.position = "bottom",
      axis.title.y = element_text(size = 6)
    ),
  filename = here::here("main/fig/fig2-comparison"),
  size = std_size$third,
  formats = c("pdf", "eps")
)

figs2 = do_figure_2(qual_df9_7_day)

.gg_save_as(
  figs2$s0 +
    figs2$main +
    plot_layout(ncol = 1, heights = c(1, 1)) +
    plot_annotation(tag_levels = "A"),
  filename = here::here("s1/fig/fig5-7-day-scenario"),
  size = std_size$full,
  formats = c("pdf", "eps")
)


figs3 = do_figure_2(qual_df9_no_lag)

.gg_save_as(
  figs3$s0 +
    figs3$main +
    plot_layout(ncol = 1, heights = c(1, 1)) +
    plot_annotation(tag_levels = "A"),
  filename = here::here("s1/fig/fig6-not-lagged-scenario"),
  size = std_size$full,
  formats = c("pdf", "eps")
)
