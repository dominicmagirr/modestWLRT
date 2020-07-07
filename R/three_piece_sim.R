three_piece_sim = function(n_c = 100,
                           n_e = 100,
                           rec_period = 12,
                           rec_power = 1,
                           rate_c_1 = log(2) / 9,
                           rate_c_2 = 0.04,
                           rate_c_3 = 0.04,
                           rate_e_1 = log(2) / 9,
                           rate_e_2 = 0.04,
                           rate_e_3 = 0.04,
                           delay_1 = 6,
                           delay_2 = 12,
                           max_cal_t  = 36,
                           n_events = NULL){
  
  if (is.null(max_cal_t) && is.null(n_events)) stop("either max_cal_t or n_events must be specified.")
  if ((!is.null(max_cal_t)) && (!is.null(n_events))) stop("one of max_cal_t and n_events must be NULL.")
  if (is.null(max_cal_t) && (n_events > n_c + n_e)) stop("number of events not reached.")
  
  # simulate recruitment times from power model:
  
  rec_c = rec_period * runif(n_c) ^ (1 / rec_power)
  rec_e = rec_period * runif(n_e) ^ (1 / rec_power)
  
  # control event times are 3-piece exponentially distributed:
  t_1_c = rexp(n_c, rate = rate_c_1)
  t_2_c = rexp(n_c, rate = rate_c_2)
  t_3_c = rexp(n_c, rate = rate_c_3)
  t_c = ifelse(t_1_c < delay_1,
               t_1_c,
               ifelse(delay_1 + t_2_c < delay_2, delay_1 + t_2_c, delay_2 + t_3_c)) 
  
  # experimental event times come from 3-piece exponential distribution:
  
  t_1_e = rexp(n_e, rate = rate_e_1)
  t_2_e = rexp(n_e, rate = rate_e_2)
  t_3_e = rexp(n_e, rate = rate_e_3)
  t_e = ifelse(t_1_e < delay_1,
               t_1_e,
               ifelse(delay_1 + t_2_e < delay_2, delay_1 + t_2_e, delay_2 + t_3_e)) 
  
  # (calendar) event times, relative to start of trial.
  
  cal_t_c = rec_c + t_c
  cal_t_e = rec_e + t_e
  
  if (is.null(max_cal_t)){
    max_cal_t <- sort(c(cal_t_c, cal_t_e))[n_events]
  }
  
  # does the patient have an event before the data cut-off:
  
  event_c = cal_t_c <= max_cal_t
  event_e = cal_t_e <= max_cal_t
  
  # if patient's event time is censored, calculate their follow-up time:
  
  obs_t_c = ifelse(event_c, t_c, max_cal_t - rec_c)
  obs_t_e = ifelse(event_e, t_e, max_cal_t - rec_e)
  
  # store in data frame with group label:
  
  df = data.frame(time = c(obs_t_c, obs_t_e),
                  event = c(event_c, event_e),
                  group = rep(c("control", "experimental"), c(n_c, n_e)))
  
  # round time to 2 dp
  df$time = round(df$time, 2)
  
  df
}
