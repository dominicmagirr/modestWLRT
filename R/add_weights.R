#' Calculate weights for a weighted log-rank test.
#' 
#' \code{add_weights} Calculate weights for a weighted log-rank test and add them to the risk table.
#' @param{risk_table} A risk table produced by the function \code{get_risk_table}.
#' @param{method} The type of weighted log-rank test.
#' \code{"fixed_c"} means that the scores \code{c} are fixed at 1 while \code{t < delay}. Thereafter, the weights \code{w} remain fixed.
#' \code{"fh"} is a Fleming-Harrington test with parameters \code{rho} and \code{gamma}.
#' \code{"step"} means that the weight \code{w} is fixed at 0 while \code{t < delay}. Thereafter, the weight is fixed at 1. 
#' @param{delay} Parameter used by the \code{"fixed_c"} and \code{"step"} methods. 
#' @param{rho} First parameter for the \code{"fh"} method.
#' @param{gamma} Second parameter for the \code{"fh"} method.
#' @return A risk table consisting of the original \code{risk_table} with 3 additional columns:
#' \code{c} The scores in the score test corresponding to events (uncensored observations).
#' \code{w} The weights for the weighted log-rank test. 
#' \code{C} The scores in the score test corresponding to censored observations.
#' @export


add_weights = function(risk_table,
                       method = "fixed_c",
                       delay = 6,
                       rho = 0,
                       gamma = 0,
                       plot_weights = FALSE){
  
  
  if (method == "fixed_c"){
    
    
    
    # Kaplan-Meier estimate of survival in the lumped data:
    risk_table$s = exp(cumsum(log(1 - risk_table$d / risk_table$n)))
    risk_table$s_minus = c(1, risk_table$s[-length(risk_table$s)])
    
    if (delay == 0) max_weight = 1
    else { 
      max_weight = 1 / max(risk_table$s[risk_table$t >= delay])
    }
    
    # Modest (delay) weights:
    risk_table$w = pmin(1 / risk_table$s_minus, max_weight) 
    risk_table$c = NA
    risk_table$C = NA
    
    ## w --> C --> c
    
    # use Condition (3) in Leton and Zuluaga:
    
    for (j in 1:length(risk_table$t)){
      
      risk_table$C[j] = -with(risk_table, sum(w[1:j] * d[1:j] / n[1:j]))
      risk_table$c[j] = risk_table$C[j] + risk_table$w[j]
      
    }
    
  }
  else if (method == "fixed_c_orig"){
  
    risk_table$c = 1
    risk_table$w = 1
    risk_table$C = 1
    
    ## c --> w --> C
    
    ## start with c = 1 up to t = delay.
    
    if (any(risk_table$t < delay)) {
      
      # find maximum event time before the end of the delay:
      max_j = max(which(risk_table$t < delay))
      
    }
    else {
      max_j = 0
    }
    
    ## Use Condition (4) in Leton and Zuluaga:
    
    if (max_j > 0){
      
      for (j in 1:max_j){
        
        if (risk_table$n[j] == risk_table$d[j]){
          
          risk_table$w[j] = 0
          
        }
        else if (j == 1){
          
          risk_table$w[j] = risk_table$n[1] / (risk_table$n[1] - risk_table$d[1])
          
        }
        else{ 
          
          risk_table$w[j] = risk_table$w[j-1] * risk_table$n[j] / (risk_table$n[j] - risk_table$d[j])
          
        }
      }
      risk_table$C[1:max_j] = risk_table$c[1:max_j] - risk_table$w[1:max_j]
    }
    ## w --> C --> c
    
    ## now keep w fixed for t > delay
    
    if (max_j < length(risk_table$t)){
      
      risk_table$w[(max_j + 1):length(risk_table$t)] = risk_table$w[max(1,max_j)]
      
      
      # use Condition (3) from Leton and Zuluaga:
      
      for (j in (max_j + 1):length(risk_table$t)){
        
        risk_table$C[j] = -with(risk_table, sum(w[1:j] * d[1:j] / n[1:j]))
        risk_table$c[j] = risk_table$C[j] + risk_table$w[j]
        
      }
      
    }
    
  }
  
  else if (method == "fh"){
    
    
    # Kaplan-Meier estimate of survival in the lumped data:
    risk_table$s = exp(cumsum(log(1 - risk_table$d / risk_table$n)))
    risk_table$s_minus = c(1, risk_table$s[-length(risk_table$s)])
    
    # Fleming-Harrington (rho, gamma) weights:
    risk_table$w = risk_table$s_minus ^ rho * (1 - risk_table$s_minus) ^ gamma
    risk_table$c = NA
    risk_table$C = NA
    
    ## w --> C --> c
    
    # use Condition (3) in Leton and Zuluaga:
    
    for (j in 1:length(risk_table$t)){
      
      risk_table$C[j] = -with(risk_table, sum(w[1:j] * d[1:j] / n[1:j]))
      risk_table$c[j] = risk_table$C[j] + risk_table$w[j]
      
    }
    
    
  }
  
  else if (method == "step"){
    
    
    # weight fixed at zero before delay, 1 after delay:
    
    risk_table$w = ifelse(risk_table$t < delay, 0, 1)
    risk_table$c = NA
    risk_table$C = NA
    
    ## w --> C --> c
    
    # use Condition (3) in Leton and Zuluaga:
    
    for (j in 1:length(risk_table$t)){
      
      risk_table$C[j] = -with(risk_table, sum(w[1:j] * d[1:j] / n[1:j]))
      risk_table$c[j] = risk_table$C[j] + risk_table$w[j]
      
    }

    
  }
  else stop("unmatched method")
  
  # round w/c/C to 2dp
  
  risk_table$w = round(risk_table$w, 2)
  risk_table$c = round(risk_table$c, 2)
  risk_table$C = round(risk_table$C, 2)
  
  
  if (plot_weights){
    
    
    ltys = c("weights" = "solid",
             "score (censored)" = "dotted",
             "score (events)" = "dashed")
    
    p = ggplot(data = risk_table) +
      geom_line(aes(x = t, y = C, linetype = "score (censored)")) +
      geom_point(aes(x = t, y = C)) +
      geom_line(aes(x = t, y = c, linetype = "score (events)")) +
      geom_point(aes(x = t, y = c)) +
      geom_line(aes(x = t, y = w, linetype = "weights")) +
      geom_point(aes(x = t, y = w)) +
      scale_linetype_manual(name = "", 
                            values = ltys, 
                            limits = c("weights",
                                       "score (events)",
                                       "score (censored)"))
    
    return(list(risk_table = risk_table,
                p = p))
    
  }
  else {
    
    risk_table
    
  }
  
}
