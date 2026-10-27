#' @title mid_run_avg
#'
#' @description Estimates the SSD with the mid-run estimates method (Levitt, 1971;
#' Wetherill & Levitt, 1965).
#'
#' @param dat a dataframe containing id, ssd, type
#' @param participant Name of the column containing the participant identifier.
#' @param type the column that reflects the current trial type (stop vs. go)
#' @param ssd Name of the column containing the current stop signal delay
#' @param stop_value prevents silent failing, coding of stop trials
#' @param every_second boolean use every second run or use all (full) runs
#' @return a dataframe containing an estimate of the SSD and the number of runs
#' used to calculate the SSD
#' @examples
#' dat <- data.frame(
#'   id = 1,
#'   trial_type = "stop",
#'   ssd = c(100, 150, 200, 150, 100, 150, 200, 150)
#' )
#' mid_run_avg(dat, participant = id, trial_type = trial_type, ssd = ssd)
#' mid_run_avg(dat, participant = id, trial_type = trial_type, ssd = ssd, every_second = FALSE)
#' @export
mid_run_avg <- function(dat,
                        participant = id,
                        trial_type = trial_type,
                        ssd = ssd,
                        stop_value = "stop",
                        every_second = TRUE) {
  return(dat |>
    dplyr::group_by({{ participant }}) |>
    dplyr::filter({{ type }} == stop_value) |>
    dplyr::filter(is.na({{ ssd }} - dplyr::lag({{ ssd }})) | {{ ssd }} != dplyr::lag({{ ssd }})) |>
    dplyr::mutate(
      step_in = {{ ssd }} - dplyr::lag({{ ssd }}),
      step_out = dplyr::lead({{ ssd }}) - {{ ssd }},
      extrema = dplyr::case_when(
        step_in > 0 & step_out < 0 ~ "max",
        step_in < 0 & step_out > 0 ~ "min",
        is.na(step_in) ~ "start",
        TRUE ~ "other"
      )
    ) |>
    dplyr::filter(extrema != "other") |>
    dplyr::mutate(
      run_id  = dplyr::row_number(),
      use_run = if (every_second) run_id %% 2 == 0 else extrema != "start",
      run_avg = dplyr::if_else(use_run, ({{ ssd }} + dplyr::lead({{ ssd }})) / 2, NA_real_)
    ) |>
    dplyr::summarise(
      mid_run_estimate = mean(run_avg, na.rm = TRUE),
      n_runs = sum(!is.na(run_avg))
    ))
}
