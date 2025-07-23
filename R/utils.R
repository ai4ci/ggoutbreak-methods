#' Format `ggoutbreak` data for `scoringutils`
#'
#' @param prediction the output from ggoutbreak `incidence` and `rt` predictions + maybe grouping
#' @param raw the raw incidence dataframe (`time` and `count`) + same grouping as prediction
#' @param changes a change point or time series for $R_t$ with columns `t` and `R_t`
#' @param ip the infectivity profile(s) with columns `boot`,`time`, and `probability`:
#'
#' @return a long format suitable for `scoringutils::score`
.scoringutils = function(prediction, raw, changes, ip) {

  omega_t = ip %>% dplyr::filter(dplyr::cur_group_id()==1) %>% dplyr::pull(probability)

  rt_true = raw %>% dplyr::transmute(
    time,
    true_value = rt,
    target_type = "rt"
  ) %>% dplyr::filter(!is.na(true_value))

  incid_true = raw %>%
    transmute(time, true_value = count, target_type="incidence")

  shared = intersect(group_vars(prediction),group_vars(raw))

  prediction %>%
    dplyr::select(c(time,tidyselect::matches(".*\\.0\\.[0-9]+"))) %>%
    tidyr::pivot_longer(-c(time,!!!groups(prediction)),names_pattern = "([^\\.]*)\\.([0-9]+\\.?[0-9]*)", names_to = c("target_type","quantile"), values_to = "prediction") %>%
    dplyr::mutate(quantile = as.numeric(quantile)) %>%
    dplyr::inner_join( dplyr::bind_rows(rt_true, incid_true), by = c("time","target_type", shared) )
}


.true_rt = function(max_time, rt_fn, ip) {

  omega_t = ggoutbreak::summarise_ip(ip) %>% dplyr::pull(probability)

  rt_true = tibble::tibble(
    time = ggoutbreak::as.time_period(1:max_time,"1 days"),
    rt_inst = rt_fn(t=1:max_time),
    rt_case = rev(as.numeric(stats::filter(rev(rt_fn(t=1:max_time)), omega_t, sides=1)))
  )

  return(rt_true)
}
