#' @title Mean estimation of the SSRT
#'
#' @description Estimates the SSRT with the mean method: mean go RT minus the
#'  mean SSD. This assumes that p(response | signal) = 0.5.
#' @importFrom rlang :=
#' @import dplyr
#' @import cli
#' @param dat a dataframe with one row per trial.
#' @param participant column containing the participant identifier.
#' @param trial_type column containing the trial type ("go" vs. "stop").
#' @param corr logical column, `TRUE` if the trial was correct.
#' @param ssd column containing the stop-signal delay.
#' @param rt column containing the reaction times.
#' @param RTexc logical; if `TRUE`, go trials with RTs outside
#'  `[RTmin, RTmax]` are excluded from the mean go RT.
#' @param RTrep logical; if `TRUE`, go go-omissions are replaced with the
#'  RTmax value.
#' @param ErrExc logical; if `TRUE`, incorrect trials are excluded from the
#'  analysis
#' @param RTmin,RTmax bounds of the RT window used when `RTexc = TRUE`.
#' @param mid_run logical; if `TRUE`, the mid-run SSD estimate is subtracted
#'  from the mean go RT, otherwise the mean SSD. The mid-run estimate is
#'  returned either way.
#' @return a dataframe with one row per participant: `meanRT`, `meanSSD`,
#'  `mid_run_estimate`, `n_runs` and `SSRT`.
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   id = 1,
#'   trial_type = rep(c("go", "stop"), c(60, 20)),
#'   corr = c(rep(TRUE, 60), rep(c(TRUE, FALSE), 10)),
#'   ssd = c(rep(NA, 60), rep(c(100, 150, 200, 150), 5)),
#'   rt = c(round(rnorm(60, 500, 80)), rep(NA, 20))
#' )
#' meanSSRT(dat,
#'   participant = id, trial_type = trial_type, corr = corr,
#'   ssd = ssd, rt = rt
#' )
#' @export
meanSSRT <- function(
  dat,
  participant = id,
  trial_type = trial_type,
  corr = corr,
  ssd = ssd,
  rt = rt,
  RTexc = FALSE,
  RTrep = FALSE,
  ErrExc = FALSE,
  RTmin = 100,
  RTmax = 600,
  mid_run = TRUE
) {
  cols <- c(
    participant = rlang::as_name(rlang::enquo(participant)),
    trial_type  = rlang::as_name(rlang::enquo(trial_type)),
    corr        = rlang::as_name(rlang::enquo(corr)),
    ssd         = rlang::as_name(rlang::enquo(ssd)),
    rt          = rlang::as_name(rlang::enquo(rt))
  )

  if (!is.data.frame(dat)) {
    cli::cli_abort("{.arg dat} must be a data frame, not {.cls {class(dat)}}.")
  }

  bad_cols <- cols[!cols %in% names(dat)]
  if (length(bad_cols) > 0) {
    rlang::abort(c(
      "Column(s) not found in `dat`.",
      x = paste0("`", names(bad_cols), "` -> `", bad_cols, "`"),
      i = paste("Available:", paste(names(dat), collapse = ", "))
    ))
  }

  by_col <- rlang::as_name(rlang::enquo(participant))
  dat_mid_run <- mid_run_avg(
    dat,
    participant = {{ participant }},
    trial_type  = {{ trial_type }},
    ssd         = {{ ssd }}
  )

  if (RTrep) {
    dat <- dat |> dplyr::mutate(
      {{ rt }} := dplyr::if_else(
        condition = is.na({{ rt }}),
        true      = RTmax,
        false     = {{ rt }}
      )
    )
  }

  if (RTexc) {
    dat <- dat |> dplyr::mutate(
      {{ rt }} := dplyr::case_when(
        {{ trial_type }} == "stop" ~ {{ rt }},
        {{ rt }} > RTmax ~ NA_integer_,
        {{ rt }} < RTmin ~ NA_integer_,
        TRUE ~ {{ rt }}
      )
    )
  }

  if (ErrExc) {
    dat <- dat |> dplyr::filter({{ trial_type }} == "stop" | {{ corr }})
  }

  dat |>
    dplyr::group_by({{ participant }}) |>
    dplyr::summarise(
      meanRT = mean(
        {{ rt }}[{{ trial_type }} == "go"],
        na.rm = TRUE
      ),
      meanSSD = mean({{ ssd }}[{{ trial_type }} == "stop"], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(dat_mid_run, by = by_col) |>
    dplyr::mutate(
      SSRT = if (mid_run) meanRT - mid_run_estimate else meanRT - meanSSD
    ) |>
    dplyr::mutate(
      meanRT  = round(meanRT),
      meanSSD = round(meanSSD),
      SSRT    = round(SSRT)
    )
}

#' @title integration method of computing the SSRT
#'
#' @description utilises the RT distribution of go trials in a
#'  Stop-Signal task to estimate the Stop-Signal Reaction Time using the
#'  distribution of go-trials (with replacement of go-omissions) and the
#'  probability of giving a response when a signal is present for each
#'  participant (see Verbruggen et al., 2019).
#' @importFrom rlang :=
#' @import dplyr
#' @import cli
#' @param dat a dataframe containing id number, rt, accuracy, the conditions
#'  and stop-signal delays
#' @param participant Name of the column containing the personal identifier.
#' @param trial_type the column that reflects the current trial type (stop vs. go)
#' @param corr boolean that reflects whether the current trial is Correct
#' @param RTmax maximum RT boundary to determine slow responses
#' @param ssd Name of the column containing the current stop signal delay
#' @param rt Reaction times
#' @param replace_slow boolean that indicates whether slow values should be
#'  replaced with the slowest correct RT, to correct distribution (see
#'  Verbruggen et al., 2019)
#' @param mid_run boolean whether the mid-run average should be used to compute
#'  the SSD
#' @return a dataframe containing the id, an estimate of the SSRT,
#'  nth reaction time, used to compute the SSRT, the mean Stop Signal delay,
#'  the p(response | signal), the mid-run estimate of the SSD, the number of
#'  runs that went into the estimation of the SSD and the estimated SSRT.
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   id = 1,
#'   trial_type = rep(c("go", "stop"), c(60, 20)),
#'   corr = c(rep(TRUE, 60), rep(c(TRUE, FALSE), 10)),
#'   ssd = c(rep(NA, 60), rep(c(100, 150, 200, 150), 5)),
#'   rt = c(round(rnorm(60, 500, 80)), rep(NA, 20))
#' )
#' intSSRT(dat,
#'   participant = id, trial_type = trial_type, corr = corr,
#'   ssd = ssd, rt = rt, mid_run = FALSE
#' )
#' @export
intSSRT <- function(
  dat,
  participant = id,
  trial_type = trial_type,
  corr = corr,
  ssd = ssd,
  rt = rt,
  replace_slow = TRUE,
  RTmax = 600,
  mid_run = TRUE
) {
  cols <- c(
    participant = rlang::as_name(rlang::enquo(participant)),
    trial_type  = rlang::as_name(rlang::enquo(trial_type)),
    corr        = rlang::as_name(rlang::enquo(corr)),
    ssd         = rlang::as_name(rlang::enquo(ssd)),
    rt          = rlang::as_name(rlang::enquo(rt))
  )

  if (!is.data.frame(dat)) {
    cli::cli_abort("{.arg dat} must be a data frame, not {.cls {class(dat)}}.")
  }

  bad_cols <- cols[!cols %in% names(dat)]
  if (length(bad_cols) > 0) {
    rlang::abort(c(
      "Column(s) not found in `dat`.",
      x = paste0("`", names(bad_cols), "` -> `", bad_cols, "`"),
      i = paste("Available:", paste(names(dat), collapse = ", "))
    ))
  }

  by_col <- rlang::as_name(rlang::enquo(participant))
  dat_mid_run <- mid_run_avg(
    dat,
    participant = {{ participant }},
    trial_type  = {{ trial_type }},
    ssd         = {{ ssd }}
  )

  if (replace_slow) {
    dat <- dat |>
      dplyr::mutate(
        slow = dplyr::if_else(
          condition = {{ rt }} > RTmax & {{ trial_type }} == "go" | is.na({{ rt }})[{{ trial_type }} == "go"],
          true      = TRUE,
          false     = FALSE
        )
      ) |>
      dplyr::group_by({{ participant }}) |>
      dplyr::mutate(
        {{ rt }} := dplyr::if_else(
          condition = slow,
          true      = max({{ rt }}[!slow & {{ trial_type }} == "go" & {{ corr }}], na.rm = TRUE),
          false     = {{ rt }}
        )
      ) |>
      dplyr::ungroup()
  }

  dat_SSRT <- dat |>
    dplyr::group_by({{ participant }}) |>
    dplyr::summarise(
      pResp   = mean(!{{ corr }}[{{ trial_type }} == "stop"], na.rm = TRUE),
      meanSSD = mean({{ ssd }}[{{ trial_type }} == "stop"], na.rm = TRUE),
      nGOs    = sum({{ trial_type }} == "go"),
      nthRT   = sort({{ rt }}[{{ trial_type }} == "go"])[max(1L, round(pResp * nGOs))],
      .groups = "drop"
    )
  dat_SSRT <-
    dplyr::left_join(
      x  = dat_SSRT,
      y  = dat_mid_run,
      by = by_col
    ) |>
    dplyr::mutate(
      SSRT = if (mid_run) nthRT - mid_run_estimate else nthRT - meanSSD
    ) |>
    dplyr::mutate(
      pResp            = round(pResp * 100, 2),
      meanSSD          = round(meanSSD),
      nGOs             = nGOs,
      nthRT            = round(nthRT),
      mid_run_estimate = round(mid_run_estimate),
      SSRT             = round(SSRT)
    )


  return(dat_SSRT)
}
