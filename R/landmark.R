#' Perform a landmark analysis at a specified (patient) time.
#' 
#' \code{landmark} performs a landmark analysis at a specified (patient) time. Based on Kaplan-Meier estimates of survival 
#' on the two treatment arms.
#' @param{df} Data frame containing simulated survival data set in standard format. 
#' Three columns: survival time \code{time}, whether patient has an \code{event} (1 = yes, 0 = censored), 
#' and treatment \code{group" (\code{control" or \code{experimental}).
#' @param{time} The (patient) times to perform the landmark analysis at. E.g., survival at \code{time = c(6,12,18)} months.
#' @return \code{z} the standardized test statistic. Large values indicate better survival on the experimental arm.
#' @export



landmark = function(df, time){

  ## use survival::survfit to do KM estimation by group:
  fit <- survival::survfit(Surv(time, event) ~ group, data = df)
  
  if(any(min(fit$time[fit$n.risk == 1]) < time)) return(rep(NA, length(time)))
  
  info = summary(fit, time = time)
  
  
  ## extract survival probabilities with standard errors:
  
  info_df = data.frame(t = info$time,
                       p = info$surv,
                       se = info$std.err,
                       group = info$strata)
  
  
  
  z_stats = info_df %>% 
              group_by(t) %>%
              summarize(diff = p[group == "group=experimental"] - p[group == "group=control"],
                        se = sqrt(sum(se ^ 2))) %>%
              mutate(z = diff/se)

  z_stats$z
  
}


