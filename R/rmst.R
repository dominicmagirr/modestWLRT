#' Perform a restricted-mean survival time.
#' 
#' \code{rmst} performs a restricted mean survival time analysis at a specified (patient) time.
#' @param{df} Data frame containing simulated survival data set in standard format. 
#' Three columns: survival time \code{time}, whether patient has an \code{event} (1 = yes, 0 = censored), 
#' and treatment \code{group" (\code{control" or \code{experimental}).
#' @param{time} The (patient) time to perform the RMST analysis at. Must be length 1.
#' @return \code{z} the standardized test statistic. Large values indicate better survival on the experimental arm.
#' @export



rmst = function(df, time){
  
  if (length(time) != 1) stop("time must be length 1")
  
  ## use survival::survfit to do KM estimation by group:
  fit <- survival::survfit(Surv(time, event) ~ group, data = df)
  
  if(any(min(fit$time[fit$n.risk == 1]) < time)) return(NA)
  
  ## use survRM2:: to do RMST estimation by group:
  
  results <- rmst2(time = df$time, 
                   status = as.numeric(df$event), 
                   arm = ifelse(df$group == "control", 0, 1),
                   tau = time)
  
  
  mean_1 <- results$RMST.arm1$rmst["Est."]
  se_1 <- results$RMST.arm1$rmst["se"]
  
  mean_0 <- results$RMST.arm0$rmst["Est."]
  se_0 <- results$RMST.arm0$rmst["se"]
  
  mean_diff <- mean_1 - mean_0
  se_diff <- sqrt(se_0 ^ 2 + se_1 ^2)
  
  z_diff <- mean_diff / se_diff
  
  z_diff
  
  
}


