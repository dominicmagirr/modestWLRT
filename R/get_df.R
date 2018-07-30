#' Simulate survival data from a two-arm trial
#'
#' \code{get_df_strong} Simulate survival data: exponential vs. two-piece exponential.
#' \param{n_c} Number of patients on control treatment.
#' \param{n_e} Number of patients on experimental treatment.
#' \param{rec_period} Recruitment period.
#' \param{rec_power} Recruitment follows a power model. Pr(recruited before T) = (T / rec_period) ^ rec_power.
#' \param{med_c} Median survival time on control.
#' \param{rate_e_1} Event rate during first period on experimental arm.
#' \param{rate_e_2} Event rate during second period on experimental arm.
#' \param{delay} Length of first period.
#' \param{max_cal_t} Maximum calendar time, i.e., time from start of the trial to data cut-off.
#' \return A data frame containing survival time, whether patient has event (1 = yes, 0 = censored), and treatment arm.
#' @export


get_df = function(n_c = 100,
                  n_e = 100,
                  rec_period = 12,
                  rec_power = 1,
                  med_c = 15,
                  rate_e_1 = log(2) / 9,
                  rate_e_2 = 0.04,
                  delay = 6,
                  max_cal_t  = 36){

  # simulate recruitment times from power model:

  rec_c = rec_period * runif(n_c) ^ (1 / rec_power)
  rec_e = rec_period * runif(n_e) ^ (1 / rec_power)

  # control event times are exponentially distributed:

  t_c = rexp(n_c, rate = log(2) / med_c)

  # experimental event times come from 2-piece exponential distribution:

  t_1_e = rexp(n_e, rate = rate_e_1)
  t_2_e = rexp(n_e, rate = rate_e_2)
  t_e = ifelse(t_1_e < delay, t_1_e, delay + t_2_e)

  # (calendar) event times, relative to start of trial.

  cal_t_c = rec_c + t_c
  cal_t_e = rec_e + t_e

  # does the patient have an event before the data cut-off:

  event_c = cal_t_c < max_cal_t
  event_e = cal_t_e < max_cal_t

  # if patient's event time is censored, calculate their follow-up time:

  obs_t_c = ifelse(event_c, t_c, max_cal_t - rec_c)
  obs_t_e = ifelse(event_e, t_e, max_cal_t - rec_e)

  # store in data frame with group label:

  df = data.frame(time = c(obs_t_c, obs_t_e),
                  event = c(event_c, event_e),
                  group = rep(c("control", "experimental"), c(n_c, n_e)))

  df
}

