#' Convert simulated survival data set into risk table format.
#' 
#' \code{get_risk_table} Calculate risk table for a simulated two-arm survival data set.
#' @param{df} Data frame containing simulated survival data set in standard format. 
#' Three columns: survival time \code{time}, whether patient has an \code{event} (1 = yes, 0 = censored), 
#' and treatment \code{group" (\code{control" or \code{experimental}).
#' @return A risk table with columns:
#' \code{t} the event times, in ascending order
#' \code{n_e} the number of patients at risk on the experimental treatment arm just prior to \code{t}.
#' \code{n_c} the number of patients at risk on the control treatment arm just prior to \code{t}.
#' \code{d_e} the number of events on the experimental arm at time \code{t}.
#' \code{d_c} the number of events on the control arm at time \code{t}.
#' \code{n} = \code{n_e} + \code{n_c}.
#' \code{d} = \code{d_e} + \code{d_c}.
#' \code{l} = \code{l_e} + \code{l_c}.
#' \code{l_e} the number of patients on the experimental treatment arm who censored after the current \code{t} but before 
#' the subsequent \code{t}.
#' @export



get_risk_table = function(df){
  
  # arrange the data set in increasing order of survival time:
  
  df = df[order(df$time),]
  
  # split into 2 data sets: one for control; one for experimental:
  
  df_c = df[df$group == "control",]
  df_e = df[df$group == "experimental",]
  
  # number of patients on each arm
  
  n_c = length(df_c$time)
  n_e = length(df_e$time)
  
  # the number of patients at risk will decrease by 1 after each event/censored observation.
  
  at_risk_c = n_c - 1:n_c + 1
  at_risk_e = n_e - 1:n_e + 1
  
  # create a risk table just for the control arm data...
  
  risk_table_c = data.frame(t = df_c$time,
                            n_c = at_risk_c,
                            d_c = as.numeric(df_c$event))
  
  # ...where there no patients/events on the experimental arm:
  
  risk_table_c$d_e = 0
  risk_table_c$n_e = NA
  
  # create a risk table just for the experimental arm data...
  
  risk_table_e = data.frame(t = df_e$time,
                            n_e = at_risk_e,
                            d_e = as.numeric(df_e$event))
  
  # ...where there are no patients/events on the control arm:
  
  risk_table_e$d_c = 0
  risk_table_e$n_c = NA
  
  # put the risk tables on top of each other...
  
  risk_table = rbind(risk_table_c, risk_table_e)
  
  # ...and reorder by event/censoring times (across both arms):
  
  risk_table = risk_table[order(risk_table$t),]
  
  # whenever is.na(n_e) == TRUE, this means that the event/censored observation on this 
  # row was from the control arm. To fill in the number at risk on the experimental
  # arm we look at the subsequent row, repeating if necessary, until we find a row
  # where is.na(n_e) == FALSE.
  # similarly for n_c when is.na(n_c) == TRUE.
  
  risk_table = risk_table %>% tidyr::fill(n_e, n_c, .direction = "up")
  
  # at the bottom of the risk table, it's still possible that is.na(n_e) == TRUE if
  # all subsequent events/censorings are from the control arm. In this case the 
  # number at risk is zero. Similarly for the control arm.
  
  risk_table$n_c[is.na(risk_table$n_c)] = 0
  risk_table$n_e[is.na(risk_table$n_e)] = 0
  
  # now we deal with ties. We group together the data that have the same value of "t",
  # work out how many patients were at risk just prior to "t", and how many events
  # happened at "t":
  
  risk_table = risk_table %>% 
    group_by(t) %>%
    summarize(n_e = max(n_e),
              d_e = sum(d_e),
              n_c = max(n_c),
              d_c = sum(d_c)) %>%
    as.data.frame()
  
  # we only keep the "t" where there was at least one event:
  
  risk_table = risk_table[risk_table$d_e > 0 | risk_table$d_c > 0,]
  
  # calculate number of events, number at risk across arms.
  
  risk_table$n = risk_table$n_e + risk_table$n_c
  risk_table$d = risk_table$d_e + risk_table$d_c
  
  # calculate the number censored between consecutive event times:
  
  risk_table$l = risk_table$n - risk_table$d - c(risk_table$n[-1], 0)
  risk_table$l_c = risk_table$n_c - risk_table$d_c - c(risk_table$n_c[-1], 0)
  risk_table$l_e = risk_table$n_e - risk_table$d_e - c(risk_table$n_e[-1], 0)
  
  # return the completed risk table:
  
  risk_table
  
}
