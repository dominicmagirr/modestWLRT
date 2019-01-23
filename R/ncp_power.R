############### model extend #######################

model_extend <- function(model,
                         max_t,
                         length_t){
  
  change_points_plus <- unique(sort(c(seq(0, 
                                          max_t, 
                                          length.out = length_t), 
                                      model$change_points)))
  
  which_lambda <- purrr::map_dbl(change_points_plus, 
                                 function(x) sum(x >= model$change_points)) + 1
  
  
  list(change_points = change_points_plus[-1],
       lambdas_0 = model$lambdas_0[which_lambda],
       lambdas_1 = model$lambdas_1[which_lambda])
  
  
}


########## \bar{S} ############

s_bar <- function(t, recruitment, model){
  
  a_0 <- recruitment$n_0 / (recruitment$n_0 + recruitment$n_1)
  a_1 <- 1 - a_0
  
  s_0 <- expectedevents:::surv_pieces_simple(t, model$change_points, model$lambdas_0)
  s_1 <- expectedevents:::surv_pieces_simple(t, model$change_points, model$lambdas_1)
  
  a_0 * s_0 + a_1 * s_1
  
}

########## w ##################

w <- function(t, t_star, recruitment, model){
  
  s_bar(pmin(t, t_star), recruitment, model) ^ (-1)
  
}


######### ncp ############

#' Find the approximate non-centrality parameter and power of a modestWLRT.
#' 
#' \code{ncp_power} returns the approximate non-centrality parameter and power of a modestWLRT for a range
#' of possible values of t*.
#' @param{t_star} A vector. A range of possible values for t*.
#' @param{model} A piecewise constant hazard model.
#'   A list containing the \code{change_points}; the rates \code{lambdas_0} on the control arm; 
#'   and the rates \code{lambdas_1} on the treatment arm.
#' @param{recruitment} List of recruitment information. 
#'   Containing \enumerate{
#'                 \item Sample size on control, \code{n_0} 
#'                 \item Sample size on treatment, \code{n_1} 
#'                 \item Recruitment period, \code{r_period}
#'                 \item Recruitment parameter for power model, \code{k} 
#'               }
#' @param{dco} Time of data cut-off. 
#' @param{length_t}  
#' @param{alpha_one_sided} 
#' @return A list containing
#' \enumerate{
#'                 \item non-centrality parameter corresponding to each t* \code{ncp} 
#'                 \item power corresponding to each t* \code{power}
#'           }
#' @export


ncp_power <- function(t_star, 
                      model,
                      recruitment,
                      dco,
                      length_t,
                      alpha_one_sided = 0.025){
  
  R <- recruitment$n_1 / recruitment$n_0
  
  model_e <- model_extend(model, max_t = dco, length_t = length_t)
  
  events_info <- expectedevents::expected_events_two_arm(dco = dco,
                                                         recruitment = recruitment,
                                                         model = model_e,
                                                         mu = 0, 
                                                         total_only = FALSE)
  
  prop_events <- events_info$prop_events
  total_events <- events_info$total_events
  
  mid_t <- c(0, model_e$change_points[-length(model_e$change_points)]) + diff(c(0, model_e$change_points)) / 2
  
  ncp <- numeric(length(t_star))
  power <- numeric(length(t_star))
  
  for (i in seq_along(ncp)){
    
    weights_t_star <- purrr::map_dbl(mid_t,
                                     w,
                                     t_star = t_star[i],
                                     recruitment = recruitment,
                                     model = model_e)
    

    weights_t_star <- c(weights_t_star, 0)
    
    log_hr <- log(model_e$lambdas_1 / model_e$lambdas_0)
    
    
    num <- sum(weights_t_star * log_hr * prop_events)
    denom <- sqrt(sum(weights_t_star ^ 2 * prop_events))
    
    ncp[i] <- num / denom * sqrt(total_events * R / (R + 1) ^ 2)
    power[i] <- pnorm(qnorm(alpha_one_sided),
                      mean = num / denom * sqrt(total_events * R / (R + 1) ^ 2))
    
  }
  list(ncp = ncp,
       power = power)
}



