#' @title meanSSRT
#'
#' @description Estimates the SSRT using the mean of all SSDs in the
#' experiment. This calculation hinges upon the assumption that p(response |
#' signal = 0.5)
#'
#' @import magrittr
#' @import dplyr
#'
#' @param dat a dataframe containing vp number, rt, accuracy, the conditions
#' and stop-signal delays
#' @param participant Name of the column containing the participant identifier.
#' @param type the column that reflects the current trial type (stop vs. go)
#' @param corr boolean that reflects wether the current trial is Correct
#' @param ssd Name of the column containing the current stop signal delay
#' @param rt Reaction times
#' @return a dataframe containing an estimate of the SSRT, the mean reaction
#' times on go trials and the mean stop signal delay for each VP
#' @examples dat <- meanSSRT(dat)
#' @export
meanSSRT <- function(
  dat,
  participant = id,
  trial_type = trial_type,
  corr = corr,
  ssd = ssd,
  rt = rt,
  RTexc = FALSE,
  RTmin = 100,
  RTmax = 600
) {
  if (RTexc) {
    dat <- dat |> dplyr::filter(between({{ rt }}, RTmin, RTmax))
  }
  return(
    merge(
      dat |>
        dplyr::filter({{ trial_type }} == "go") |>
        group_by({{ participant }}) |>
        summarise(meanRT = round(mean({{ rt }}))),
      dat |>
        dplyr::filter({{ trial_type }} == "stop") |>
        group_by() |>
        summarise(meanSSD = round(mean({{ ssd }}))),
      by = {{ participant }}
    ) |> mutate(SSRT = meanRT - meanSSD)
  )
}

#' @title intSSRT
#'
#' @description utilises the RT distribution of correct Go Trials in a
#' Stop-Signal task to estimate the Stop-Signal Reaction Time
#'
#' @import magrittr
#' @import dplyr
#'
#' @param dat a dataframe containing vp number, rt, accuracy, the conditions
#' and stop-signal delays
#' @param vp_num Name of the column containing the vp identifier.
#' @param type the column that reflects the current trial type (stop vs. go)
#' @param corr boolean that reflects wether the current trial is Correct
#' @param ssd Name of the column containing the current stop signal delay
#' @param rt Reaction times
#' @param replace_slow boolean that indicates whether slow values should be
#' raplaced with the slowest correct RT, to correct distribution (see
#' Verbruggen et al., 2019)
#' @return a dataframe containing the vp identifier, an estimate of the SSRT,
#' nth reaction time, used to compute the SSRT, the mean Stop Signal delay,
#' the p(response | signal) and the estimated SSRT.
#' @examples dat <- intSSRT(dat)
#' @export
intSSRT <- function(
  dat,
  participant = id,
  trial_type = trial_type,
  corr = corr,
  ssd = ssd,
  rt = rt,
  slow = slow,
  replace_slow = TRUE
) {
  if (replace_slow) {
    dat <- dat |>
      dplyr::group_by({{ participant }}) |>
      dplyr::mutate(
        {{ rt }} := ifelse(
          test = {{ slow }},
          yes = max({{ rt }}[!{{ slow }} & {{ trial_type }} == "go"], na.rm = TRUE),
          no = {{ rt }}
        )
      ) |>
      dplyr::ungroup()
  }

  dat_SSRT <- dat |>
    dplyr::group_by({{ participant }}) |>
    dplyr::summarise(
      pResp   = mean(!{{ corr }}[{{ trial_type }} == "stop"]),
      meanSSD = mean({{ ssd }}[{{ trial_type }} == "stop"], na.rm = TRUE),
      nGOs    = sum({{ trial_type }} == "go"),
      nthRT   = sort({{ rt }}[{{ trial_type }} == "go"])[max(1L, as.integer(pResp * nGOs))],
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pResp   = round(pResp * 100, 2),
      meanSSD = round(meanSSD),
      nGOs    = nGOs,
      nthRT   = round(nthRT),
      SSRT    = nthRT - meanSSD
    )
  return(dat_SSRT)
}
